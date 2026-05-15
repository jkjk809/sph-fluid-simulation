# sph-fluid-simulation
SPH fluid sim developed with c++ and SFML. Click to watch :)

First clone Repo then cd into it.
exe file was compiled on Mac (homebrew). requires SFML
- brew install sfml@2 libomp  


 
Build Command (MacOS)

- g++ -std=c++17 sph_fluid.cpp -o sph_fluid \
-I"$(brew --prefix sfml@2)/include" \
-I"$(brew --prefix libomp)/include" \
-L"$(brew --prefix sfml@2)/lib" \
-L"$(brew --prefix libomp)/lib" \
-Xpreprocessor -fopenmp -lomp \
-lsfml-graphics -lsfml-window -lsfml-system

Run
- ./sph_fluid

Requires SFML 2.x installed via Homebrew.
Uses OpenMP for multithreaded physics acceleration.
Built with C++17, SFML2, OpenMP. Tested on macOS.

[![Watch Demo](https://img.youtube.com/vi/tVH68QJgglI/maxresdefault.jpg)](https://www.youtube.com/watch?v=tVH68QJgglI)


