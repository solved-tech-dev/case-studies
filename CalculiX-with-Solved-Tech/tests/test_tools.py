# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 解得科技
"""Local preparation/parser checks; these do not claim a solver run."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from input_deck import parse
from prepare_small import prepare

class PreparationTests(unittest.TestCase):
    def test_all_five_original_inputs(self):
        expected={'dam':(1045,320),'ax6ht':(85,32),'capacitor':(108,10),'oneel20fi':(20,1),'pan':(3745,1730)}
        with tempfile.TemporaryDirectory() as folder:
            for case,counts in expected.items():
                path=Path(folder)/(case+'.inp');prepare(case,path);model=parse(path)
                self.assertEqual((len(model['nodes']),len(model['elements'])),counts)
                self.assertEqual(model['nsets']['ALLN'],set(model['nodes']))
                self.assertEqual(model['esets']['ALLE'],set(model['elements']))
                self.assertIn('*HEAT TRANSFER',path.read_text())
                for line in path.read_text().splitlines():
                    if line and all(c in '0123456789,' for c in line):self.assertLessEqual(len(line.split(',')),16)
                with self.assertRaises(FileExistsError):prepare(case,path)

    def test_fixed_inputs_and_provenance(self):
        for case in ('dam','ax6ht','capacitor','oneel20fi','pan'):
            params=json.loads((ROOT/'cases'/case/'parameters.json').read_text())
            self.assertTrue((ROOT/'cases'/case/params['input']).is_file())
            self.assertGreater(params['measured_system']['n'],300000)
            self.assertEqual(params['solver']['ordering'],'L0')

if __name__=='__main__':unittest.main()
