# AGENTS.md

This repository is primarily developed in `python/` and non-submodule C++ source under `src/`. Keep changes focused, prefer the existing project workflow, and ask before changing build, packaging, dependency, or environment configuration.

## Python Environment

- Use the `rxgaming` conda environment for Python work.
- Do not assume the current shell or app session has the correct environment activated.
- For direct Python execution, use `conda run -n rxgaming python ...`.
- If you intentionally use plain `python`, first verify the active interpreter is from the `rxgaming` environment.
- Match the existing VS Code debug environment when relevant:
  - working directory: repo root
  - `PYTHONPATH=python`
  - `PYTHONDONTWRITEBYTECODE=1`

## Validation

- Prefer the existing Python unit tests over ad hoc smoke tests when validating Python changes.
- Python-only changes: run the relevant unit tests in `python/tests` when practical.
- C++ changes: run `cmake --build build/release --config Release`.
- Mixed Python and C++ changes: do both when practical.
- If validation is skipped or partial, say so clearly and explain why.

## Allowed Edits

- Source changes are primarily expected in `python/` and non-submodule code under `src/`.
- The following additional files are allowed to edit when directly relevant:
  - `.vscode/`
  - `.gitignore`
  - `.gitmodules`
  - `CMakeLists.txt`
  - `CMakePresets.json`
  - `environment.yml`
  - `README.md`
  - `rxgaming.spec`
  - `vcpkg.json`
  - `TODO.md`
  - `AGENTS.md`

## Ask First

- Ask before modifying submodules or proposing submodule updates.
- Ask before changing any build environment, packaging, dependency, or toolchain configuration.
- Ask before editing:
  - `CMakeLists.txt`
  - `CMakePresets.json`
  - `environment.yml`
  - `rxgaming.spec`
  - `vcpkg.json`
  - `.gitmodules`
  - `.gitignore`
- When asking, justify why the change is needed.

## Do Not Change Without Discussion

- Do not change `vcpkg` settings without prior discussion.
- Do not upgrade conda dependencies without prior discussion.
- Do not edit large resource files unless explicitly asked.
- Do not change packaging unless explicitly asked.
- Do not change the build environment without a strong reason and prior discussion.

## Submodules

- Treat these paths as submodules and do not edit them directly without approval:
  - `src/lapisgis`
  - `src/lico`
  - `src/processedfolder`
  - `src/rxtools`
- If work appears to require a submodule change, explain the needed change and ask the user to handle it or approve it first.

## Project Workflow Notes

- The standard C++ build task is `cmake --build build/release --config Release`.
- The main Python entrypoint is `python/__main__.py`.
- The repo README expects VS Code to be opened from an activated `rxgaming` prompt, but agents should still prefer `conda run -n rxgaming` for reliability.
