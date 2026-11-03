# Upstream sources

[English](README.md) | [简体中文](README.zh-CN.md) · [Home](../README.md)

`upstream/ccx_2.23.src.tar.bz2` contains the complete CalculiX 2.23 source distribution. The build tool extracts and compiles this archive. Selected source files under `upstream/src/` provide context for the solver interface and optional timing patch.

Source: [CalculiX downloads](https://www.dhondt.de/). Copyright, authorship, and GPL notices are included with the source.

Example inputs and include files are stored under `cases/<case>/source/`. Each case's `parameters.json` records its source URL, revision, mesh parameters, and measured system dimensions.

PanguLU and CUDA are external dependencies; their SDKs and licenses are installed separately.
