# Dimensions of the box

Depth/width/height offsets for a certain panel are defined as the distance
between the inner surface to the edge of the panel in the given direction.
A positive number means that the panel is shorter than the distance between
the inner surfaces.

## General

|Symbol|Value / mm|Description|
|------|----------|-----------|
|     A|    132.59|Outer height of the front, back and side panels|
|     B|    424.94|Outer width of the box|
|     C|    254.00|Outer depth of the side panels|
|     F|      2.03|Panel thickness|

## Top/bottom panels

|Symbol|Value / mm|Description|
|------|----------|-----------|
|     J|    417.83|Width (B - 7.11)|
|      |     0.375|Width offset (7.11 / 2 - 3.18)|
|     K|    238.00|Depth (C - 16.00)|
|      |     0.375|Depth offset (16.00 / 2 - 2.54 - F)|

## Front/back panels

|Symbol|Value / mm|Description|
|------|----------|-----------|
|     D|    418.84|Width (B - 6.10)|
|      |     -0.13|Width offset (6.10 / 2 - 3.18); (should be 0)|
|     A|    132.59|Height|
|      |     -1.91|Height offset|

## Side panels

|Symbol|Value / mm|Description|
|------|----------|-----------|
|     C|    254.00|Depth|
|      |     -4.57|Depth offset (-2.54 - F)|
|     A|    132.59|Height|
|      |     -1.91|Height offset|

## References

[Gray box B style data](doc/graybox_b_techdata.pdf).
We are using a full width, 3U, 10inches depth box.
