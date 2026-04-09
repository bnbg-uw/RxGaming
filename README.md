RxGaming Tool v1.0.3
Manual is available here:
https://docs.google.com/document/d/1Y\_nr5oXO7UbGjxRbcr00cXnC7lWk71\_8q\_Ht6XfXHIE/edit?tab=t.0



\## Development Setup



1\. Install prerequisites: VS Code, Anaconda, CMake, Visual Studio (C++ tools), vcpkg

2\. Set environment variable: `VCPKG\_ROOT=C:\\vcpkg`

3\. Clone repo and init submodules:

&#x20;  git clone <repo>

&#x20;  cd RxGaming

&#x20;  git submodule update --init --recursive

4\. Create conda env: `conda env create -f environment.yml`

5\. Download FIA csv's for relevant testing areas and place inside resources\\fia.

5\. Always open VS Code from an activated conda prompt:

&#x20;  conda activate rxgaming

&#x20;  code .

6\. In VS Code: CMake: Select Configure Preset -> release, then CMake: Configure

7\. The default VS Code build/debug tasks assume a Release native build of `rxgaming_core`.

8\. For `proj.db`, this project uses the vcpkg copy that matches the linked C++ PROJ library.
   CMake will auto-detect it from `VCPKG_ROOT` and `VCPKG_TARGET_TRIPLET`, or you can override
   with `-DPROJ_DB_PATH=...`.

9\. A starter PyInstaller spec is provided at `rxgaming.spec`. It explicitly bundles:
   - `python/rxgaming_core*.pyd`
   - native DLLs from `build/release/Release`
   - native DLLs from `VCPKG_ROOT/installed/<triplet>/bin`
   - the `resources/` folder

