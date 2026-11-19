"""Read-only validation of these native FEBio 4.13 uncompressed v3.5 plots.

Format follows installed FEBioPlotFile.h/.cpp, PltArchive.h and fecore_enum.h.
Reject unsupported formats instead of inferring their layout. Match converged
physical times and native node/domain/element identities, never Newton indices.
"""

import argparse
import json
import struct
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np


def chunks(data):
    offset = 0
    while offset < len(data):
        assert offset + 8 <= len(data), "Truncated chunk header"
        tag, length = struct.unpack_from("<II", data, offset)
        end = offset + 8 + length
        assert end <= len(data), f"Truncated chunk {tag:x}"
        yield tag, data[offset + 8 : end]
        offset = end
    assert offset == len(data)


def unique(data):
    result = {}
    for tag, value in chunks(data):
        assert tag not in result, f"Duplicate singleton chunk {tag:x}"
        result[tag] = value
    return result


def integer(data):
    assert len(data) == 4
    return struct.unpack("<i", data)[0]


def name(data):
    # Dictionary names are fixed-size C buffers; bytes after NUL are uninitialized.
    return bytes(data).split(b"\0", 1)[0].decode("utf-8")


def read_plot(path, input_path):
    source = ET.parse(input_path).getroot()
    control = source.find("Control")
    assert control is not None
    end = float(control.findtext("time_steps")) * float(control.findtext("step_size"))
    expected_final_time = end
    expected_nodes = sum(len(block) for block in source.findall("Mesh/Nodes"))
    expected_elements = sum(len(block) for block in source.findall("Mesh/Elements"))
    requested = {var.get("type") for var in source.findall("Output/plotfile/var")}
    del source
    raw = memoryview(Path(path).read_bytes())
    assert bytes(raw[:4]) == b"BEF\0", "Unexpected magic/endianness"
    top = list(chunks(raw[4:]))
    assert [x[0] for x in top[:2]] == [0x01000000, 0x01040000]
    root = unique(top[0][1])
    header = unique(root[0x01010000])
    assert integer(header[0x01010001]) == 0x35, "Only native version 3.5 supported"
    assert (
        integer(header[0x01010004]) == 0
    ), "Compressed plot requires a different reader"
    dictionary = {}
    regions = {
        0x01021000: 0x02020100,
        0x01023000: 0x02020300,
        0x01024000: 0x02020400,
        0x01025000: 0x02020500,
        0x01026000: 0x02020600,
    }
    for region, data in chunks(root[0x01020000]):
        assert region in regions
        for index, (tag, payload) in enumerate(chunks(data), 1):
            assert tag == 0x01020001
            item = unique(payload)
            kind, fmt = integer(item[0x01020002]), integer(item[0x01020003])
            size = integer(item[0x01020005])
            components = {0: 1, 1: 3, 2: 6, 3: 3, 4: 21, 5: 9, 6: size, 7: 3 * size}[
                kind
            ]
            assert components > 0
            dictionary[(regions[region], index)] = (
                name(item[0x01020004]),
                kind,
                fmt,
                size,
                components,
            )
    mesh = unique(top[1][1])
    nodes = unique(mesh[0x01041000])
    node_header = unique(nodes[0x01041100])
    count = integer(node_header[0x01041101])
    assert integer(node_header[0x01041102]) == 3
    node_data = nodes[0x01041200]
    assert len(node_data) == 16 * count
    ids = np.frombuffer(node_data, dtype="<i4").reshape(count, 4)[:, 0]
    coordinates = np.frombuffer(node_data, dtype="<f4").reshape(count, 4)[:, 1:]
    assert (
        np.array_equal(ids, np.arange(1, count + 1)) and np.isfinite(coordinates).all()
    )
    assert count == expected_nodes
    domains = {}
    for index, (tag, payload) in enumerate(chunks(mesh[0x01042000]), 1):
        assert tag == 0x01042100
        dom = unique(payload)
        dom_header = unique(dom[0x01042101])
        element_type = integer(dom_header[0x01042102])
        part_id = integer(dom_header[0x01042103])
        ne = integer(dom_header[0x01032104])
        connectivity = []
        for tag, element in chunks(dom[0x01042200]):
            assert tag == 0x01042201 and len(element) % 4 == 0
            connectivity.append(np.frombuffer(element, dtype="<i4"))
        assert len(connectivity) == ne and ne > 0
        elements = np.stack(connectivity)
        assert np.unique(elements[:, 0]).size == ne
        assert np.all((elements[:, 1:] >= 0) & (elements[:, 1:] < count))
        domains[index] = dict(
            type=element_type,
            part=part_id,
            elements=elements,
            nodes=np.unique(elements[:, 1:]).size,
        )
    elements_total = sum(d["elements"].shape[0] for d in domains.values())
    assert elements_total == expected_elements
    assert requested <= {
        value[0] for value in dictionary.values()
    }, "Missing requested output fields"
    states, all_headers = {}, []
    for tag, payload in top[2:]:
        assert tag == 0x02000000, f"Unsupported top-level section {tag:x}"
        state = unique(payload)
        head = unique(state[0x02010000])
        assert len(head[0x02010002]) == 4
        time = struct.unpack("<f", head[0x02010002])[0]
        status = integer(head[0x02010003])
        assert np.isfinite(time)
        fields = {}
        for region, region_data in chunks(state[0x02020000]):
            for tag, var_data in chunks(region_data):
                assert tag == 0x02020001
                var = unique(var_data)
                vid = integer(var[0x02020002])
                key = (region, vid)
                assert key in dictionary
                _, _, fmt, _, components = dictionary[key]
                for domain_id, values in chunks(var[0x02020003]):
                    values = np.frombuffer(values, dtype="<f4").astype(np.float64)
                    assert np.isfinite(values).all()
                    if region == 0x02020300:
                        assert domain_id == 0 and values.size == count * components
                    elif region == 0x02020400:
                        d = domains[domain_id]
                        ne, node_per_element = (
                            d["elements"].shape[0],
                            d["elements"].shape[1] - 1,
                        )
                        items = {0: d["nodes"], 1: ne, 2: ne * node_per_element, 3: 1}[
                            fmt
                        ]
                        assert values.size == items * components
                    elif region == 0x02020100:
                        assert domain_id == 0 and values.size == components
                    else:
                        raise AssertionError("Unsupported surface or edge output")
                    field_key = (*key, domain_id)
                    assert field_key not in fields
                    fields[field_key] = values.reshape(-1, components)
        missing = set(dictionary) - set((k[0], k[1]) for k in fields)
        assert (
            not missing
        ), f"Unexpected missing native variables at time{time}: {[dictionary[k][0] for k in sorted(missing)]}"
        flags = unique(state[0x02030000])[0x02030001]
        assert len(flags) == 4 * elements_total
        all_headers.append(
            dict(
                time=time,
                status=status,
                value_count=sum(a.size for a in fields.values()),
            )
        )
        if status == 0:
            assert time not in states, "Ambiguous repeated converged physical time"
            states[time] = dict(fields=fields, mesh_flags=bytes(flags))
    # Respect the original native output policy; pure-fluid originals explicitly
    # set plot_zero_state=0. Never turn it on or synthesize an initial field.
    if control.findtext("plot_zero_state", "1").strip() != "0":
        assert 0.0 in states
    assert states and np.isclose(max(states), expected_final_time, rtol=1e-6, atol=1e-7)
    assert (
        np.isclose(all_headers[-1]["time"], expected_final_time, rtol=1e-6, atol=1e-7)
        and all_headers[-1]["status"] == 0
    )
    return dict(
        nodes=count,
        elements=elements_total,
        ids=ids,
        coordinates=coordinates,
        domains=domains,
        dictionary=dictionary,
        states=states,
        state_headers=all_headers,
    )


def compare_plots(reference, candidate, tolerance=1e-6):
    assert (
        reference["dictionary"] == candidate["dictionary"]
    ), "Different field metadata"
    assert np.array_equal(reference["ids"], candidate["ids"])
    assert np.array_equal(reference["coordinates"], candidate["coordinates"])
    assert reference["domains"].keys() == candidate["domains"].keys()
    for index, a in reference["domains"].items():
        b = candidate["domains"][index]
        assert a["type"] == b["type"] and a["part"] == b["part"]
        assert np.array_equal(
            a["elements"], b["elements"]
        ), "Different native element identities/topology"
    assert (
        reference["states"].keys() == candidate["states"].keys()
    ), "Different converged physical-time coverage"
    results = []
    for time, a in reference["states"].items():
        b = candidate["states"][time]
        assert a["mesh_flags"] == b["mesh_flags"], "Different active element state"
        assert a["fields"].keys() == b["fields"].keys()
        for key, values in a["fields"].items():
            other = b["fields"][key]
            assert values.shape == other.shape
            difference = float(
                np.max(np.abs(other - values) / np.maximum(1.0, np.abs(values)))
            )
            results.append(
                dict(
                    time=time,
                    field=reference["dictionary"][key[:2]][0],
                    region=hex(key[0]),
                    domain_id=key[2],
                    components=values.shape[1],
                    values=int(values.size),
                    max_scaled_difference=difference,
                    passed=difference <= tolerance,
                )
            )
    populated = set((x["region"], x["field"]) for x in results)
    omitted = [
        item[0]
        for key, item in reference["dictionary"].items()
        if (hex(key[0]), item[0]) not in populated
    ]
    return dict(
        passed=all(x["passed"] for x in results),
        tolerance=tolerance,
        dictionary_fields=len(reference["dictionary"]),
        nodes=reference["nodes"],
        populated_fields=len(populated),
        native_unpopulated_variables_with_no_contact_or_rigid_entities=omitted,
        elements=reference["elements"],
        physical_times=sorted(reference["states"]),
        decoded_values_compared=sum(x["values"] for x in results),
        maximum_scaled_difference=max(x["max_scaled_difference"] for x in results),
        fields=results,
        correspondence="native node/element/domain identities, exact reference mesh, converged physical time; no Newton ordinal matching",
        precision="native xplt float32",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path, help="FEBio model used for both runs")
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    args = parser.parse_args()
    if not __debug__:
        parser.error("Run without Python optimization to enable validation")
    result = compare_plots(
        read_plot(args.reference, args.input), read_plot(args.candidate, args.input)
    )
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["passed"] else 1)
