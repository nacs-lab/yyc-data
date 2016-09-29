# The power supply used for this version of the FPGA box

## Voltages needed

1. 12V (~3A)

    Cooling fans, DAC output amplifier, DDS output amplifier

2. +-5V

    DAC output, DAC digital circuit

3. 3.3V, 1.8V

    DDS Analog and digital circuit

## Sources

The FPGA board provides a 12V power which can output no more than 5A of current.
This is shared with the FPGA board which might consume ~1A or more current
(unclear from datasheets) so it is likely not enough for all the 12V current
we need in all conditions.

The additional power for the old version comes from a 100W 48V benchtop
power supply. Preferably we should use the same input so the parts are more
usuable. The final power supplies used are,

1. 110VAC -> 48VDC

    [DTM165PW480C](http://www.digikey.com/product-detail/en/tdk-lambda-americas-inc/DTM165PW480C/285-2408-ND/5304305) 165W
    The original one was discontinued and this is the replacement
    with high power.

2. 48VDC -> +-12VDC

    [PTB48510BAS](http://www.digikey.com/product-detail/en/texas-instruments/PTB48510BAS/296-32710-ND/1573604) 72W

3. 48VDC -> 1.8VDC, 3.3VDC

    [QD48T018033-PAA0G](http://www.questcomp.com/questdetails.aspx?pn=QD48T018033-PAA0G&pnid=510949&pt=0) 77W
    This part is discontinued according to some sites.
    Bought 7 so that we can keep using this for some time.

4. 12VDC -> +-5VDC

    [PYB20-Q24-D5-U](http://www.digikey.com/product-detail/en/cui-inc/PYB20-Q24-D5-U/102-3243-ND/4477501) 20W
