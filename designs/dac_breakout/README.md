# Cable type

When connecting the daughter boards to the breakout board, there are multiple signal and power
lines that we need to connect. For each output board, we need at least,

1. +12V
2. -12V
3. +5V
4. -5V
5. GND
6. CLK
7. CS
8. DATA

(This already combines digital and analog powers.)

A few cable options,

1. Ethernet cable

    This have exactly 8 wires so it's just enough for one board per cable.
    The clock and the data should be wired in pair with a power or ground to reduce noise.
    This allows very long cables though it is questionable if the integrity of the signal
    is still good after such a long distance.

2. USB Type C cable

    This have 10x2 wires other than ground though only 10 are distinguishable naively.
    The main disadvantage is that the cable can only be 1m long max although this can possibly
    support more than one board per cable.q

    * Without direction detection

        1. A1+A12+B1+B12: GND
        2. A2+B2: CLK
        3. A3+B3: +5V
        4. A4+B4: +12V
        5. A5+B5: -12V
        6. A6+B6: CS1
        7. A7+B7: CS2
        8. A8+B8: -12V
        9. A9+B9: +12V
        10. A10+B10: -5V
        11. A11+B11: DATA

    * With direction detection

        1. A1+A12+B1+B12: GND
        2. A2: CLK
        3. B2: DATA
        4. A3+B3: +5V
        5. A4+B4: +12V
        6. A5/B5: direction detection
        7. A6+A7+B6+B7: -12V
        8. A8/B8: direction detection
        9. A9+B9: +12V
        10. A10+B10: -5V
        11. A11: CS1
        12. B11: CS2
