## inverse matrix recursive and gauss jordan methods

The whole idea of this code is to show how fast each method does the work. As you can see after running the code, gauss jordan method does it way faster than recursive method.

Using the code uploaded, with my system, the times are:
- 33993837 microseconds for recursive
- 1901 microseconds for gauss jordan

The matrix dimensions can be changed.

To run the code:
1. First, build a folder named `build` using `mkdir build`
2. Then, inside the `build` folder, run `cmake ..` in terminal
3. After that, run `make`
4. And at the end, run `./NAME-OF-THE-PROJECT`