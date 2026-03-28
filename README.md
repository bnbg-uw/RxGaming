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

5\. Always open VS Code from an activated conda prompt:

&#x20;  conda activate rxgaming

&#x20;  code .

6\. In VS Code: CMake: Select Configure Preset → debug, then CMake: Configure

