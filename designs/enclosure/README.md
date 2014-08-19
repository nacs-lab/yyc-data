# Design files for enclosure, including DDS

## Files

1. `*.fpd`

    Front Panel Designer files.

2. `*.dxf`

    For the waterjet.

3. `*.svg` and `*.png`

    Drawing

4. `dimensions`

    Dimensions of the front/back, top/bottom panels and their offsets from
    the inner corner.

## Holes and Considerations

1. Front

    1. Ethernet and USB serial console input
    2. 44V power input
    3. Mount holes for the DC-DC converter
    4. Board power
    5. SD card and switch. Leave enough space for fingers
    6. Front fan and holes next to the fan for better air flow
    7. Auxiliary USB power input

2. Back

    1. Backup ethernet and USB serial console input
    2. Backup auxiliary USB power input
    3. Back fan and holes next to the fan for better air flow

3. Top

    1. Mount holes for power splitters (8 ports x3, 3 ports x1) with different
    orientation and position options
    2. TTL/camera link with mount holes for the camera link adapter plate.
    3. CLock output
    4. REF CLK inputs and signal outputs for the DDS boards. Output holes
    should not make contact with the connectors
    5. Mount holes for the L-brackets for the DDS boards
    6. Auxiliary hole for powering the clock amplifier outside.

4. Bottom

    1. Mount holes for the FPGA board and the breakout boards.

## Drawing

1. Front

    ![Front panel](front.png)

2. Back

    ![Back panel](back.png)

3. Top

    ![Top panel](top.png)

4. Bottom

    ![Bottom panel](bottom.png)

## Parts

1. Ethernet feedthrough

    [EHRJ45P5ES <img src="img/EHRJ45P5ES.jpg" alt="EHRJ45P5ES Image" style="width: 100px;"/>](http://www.digikey.com/product-detail/en/EHRJ45P5ES/EHRJ45P5ES-ND/2666475)

2. USB feedthrough

    [17-200001 ![17-200001 Image](img/17-200001.jpg)](http://www.digikey.com/product-detail/en/17-200001/626-1352-ND/2184932)

3. L bracket

    [621 ![621 Image](img/621.jpg)](http://www.digikey.com/product-detail/en/621/621K-ND/316544)

4. Finger guard

    [08174 ![08174 Image](img/08174.jpg)](http://www.digikey.com/product-detail/en/08174/CR220-ND/43240)

5. Fan

    [HA92251V4-000U-999 <img src="img/HA92251V4-000U-999.jpg" alt="HA92251V4-000U-999 Image" style="width: 100px;"/>](http://www.digikey.com/product-detail/en/HA92251V4-000U-999/259-1614-ND/1937331)
