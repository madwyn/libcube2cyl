Cube2Cyl
========

A panoramic lib, cubic to cylindrical projection conversion

Read more details of implementation here: http://www.wenyanan.com/index.php/cube2cyl/

How to use
----------

The Cube2Cyl library is a single header file library, written in C++. There is a demo application included as usage example.

The library takes the following parameters:

    Width of the panorama (must be lager than 1)
    Height of the panorama (must be lager than 1)
    Width of the square input (the input must have a ratio of 1:1)
    View portion horizontally, [0.01, 2π]
    View portion vertically, [0.01, π]

This library takes any valid parameter and deals the ratio changes automatically. You can use super sampling to generate accurate panorama.

For the maximum performance, using a precomputed map is recommended rather than doing the calculation for each frame.
