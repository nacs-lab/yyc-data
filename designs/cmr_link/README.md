# Orientations

Both connectors (Camera link an the 2x30 connector to the breakout board) and the
decoder chip are on the same (front) side of the board, which is the side that
is closer to the FMC connector and the DDS power supply connectors on the
breakout board. Also the side that is farther from the FPGA board/chip.

When looking from the front side of the board with the 2x30 connector on the
bottom and the camera link connector on the top, the PIN 1 of the 2x30 connector
is on the bottom left of the board and it should be connected to the ground on
the breakout board as well as all other odd number pins on this connector.
Pin 60 of the connector is VCC3V3 from the FPGA board.

# Camera link pin assignments

According to [Andor's document](http://www.andor.com/learning-academy/camera-link-output-ixon-ultra-output-for-direct-data-access), each pixel is a 16bit little endian integer, and according to the standard Port A and Port B should be used in this case.

# Microstrip

Trace width: 10 mil

Trace separation: 6 mil

Trace thickness: 1.7 mil

Dielectric thickness: 11.9 mil

Relative dielectric constant: 4.5

Differential impedance: 100.4 Ohm

# Useful documents

* Sunstone

    [Sunstone capabilities](http://www.sunstone.com/pcb-capabilities/pcb-manufacturing-capabilities)

    [Sunstone board dimensions](http://www.sunstone.com/pcb-products/pcbexpress-quickturn/quickturn-pcb-construction)

* Calculators

    [Microstrip calculator](http://www.multek.se/engelska/engineering/pcb-structures-2/differential-microstrip-impedance-calculator-2)
