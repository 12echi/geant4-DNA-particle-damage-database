# Geant4-DNA damage-clustering prototype

This repository contains a Geant4-DNA simulation application with physical and chemical-stage tracking plus clustering of candidate DNA strand-break sites.

Despite the repository's historical name, the current tree is source code and example output, not a curated particle-damage database or a validated reference dataset.

## Provenance and scope

The application under `B1/` is derived from Geant4 examples. Source headers identify the Geant4-DNA collaboration and the DNA clustering example, including a DBSCAN-based clustering implementation. Preserve those notices and follow the upstream citation and licensing requirements.

The repository history does not document which upstream Geant4 release was used or which files were locally modified. Results should not be treated as reproducible until that provenance and the exact Geant4 version are recorded.

## Layout

- `B1/exampleB1.cc`: application entry point.
- `B1/src/`, `B1/include/`: detector, chemistry, tracking, and damage-clustering implementation.
- `B1/*.mac`, `B1/*.in`: example run and visualization macros.
- `B1/exampleB1.out`: example output only; it is not validation evidence.
- `B1/README`: inherited Geant4 B1 instructions, some of which may not describe the modified application accurately.

## Build

Install a Geant4 build that includes the required UI/visualization components and Geant4-DNA physics and chemistry data. From a separate build directory:

```bash
cmake -S B1 -B build
cmake --build build -j
```

For a batch-only build:

```bash
cmake -S B1 -B build -DWITH_GEANT4_UIVIS=OFF
cmake --build build -j
```

## Run

The CMake build copies the macro and input files into the build directory. Example invocations are:

```bash
cd build
./exampleB1 run2.mac
./exampleB1 exampleB1.in
```

Parameter defaults and damage classification logic are implemented in the source. Record the Geant4 version, physics and chemistry configuration, random seeds, macro files, clustering parameters, and normalization before comparing results.

## Known gaps

- No automated regression or reference-result comparison is included.
- Upstream version and local modification history are undocumented.
- The repository does not contain a structured damage database.
- `B1/test.cpp` is empty.
- The checked-in example output cannot validate a rebuilt executable.

## License

No repository-level license has been declared. The source contains upstream Geant4-DNA notices; do not add a conflicting license or assume redistribution rights until provenance is established.
