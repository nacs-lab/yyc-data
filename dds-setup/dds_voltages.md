|Position easured                 |Value / V |
|---------------------------------|----------|
|DC-DC converter in, on board     |47.98     |
|3.3V out, closest point to converter out|3.610|
|3.3V out, furthest point to converter out on board|3.602|
|1.8V out, closest point to converter out|1.966|
|1.8V out, furthest point to converter out on board|1.957|
|Voltage drop for 1.8V line on board|0.004   |
|Voltage drop for 3.3V line on board|0.003   |
|Voltage drop for converter ground out on board|0.005|
|3.3V at the connector on the breakout board 1|3.464|
|3.3V at the connector on the breakout board 2|3.412|
|1.8V at the connector on the breakout board 1|1.832|
|1.8V at the connector on the breakout board 2|1.817|
|3.3V at the jumper on the breakout board 1|3.417|
|3.3V at the jumper on the breakout board 2|3.384|
|1.8V at the jumper on the breakout board 1|1.809|
|1.8V at the jumper on the breakout board 2|1.803|
|3.3V at the 8x2 connector on the DDS board 1|3.423 (expect [3.135, 3.465])|
|3.3V at the 8x2 connector on the DDS board 2|3.390 (expect [3.135, 3.465])|
|1.8V at X6 on the DDS board 1|1.800 (expect [1.71, 1.89])|
|1.8V at X6 on the DDS board 2|1.796 (expect [1.71, 1.89])|
|3.3V total voltage drop for board 1|0.187   |
|3.3V total voltage drop for board 2|0.226   |
|1.8V total voltage drop for board 1|0.166   |
|1.8V total voltage drop for board 2|0.169   |
|12V on board                     |11.96     |
|12V on breakout board 1          |11.84     |
|12V on breakout board 2          |11.86     |
|Output amp bias on breakout 1|4.446 (expect [4.3, 4.9])|
|Output amp bias on breakout 2|4.413 (expect [4.3, 4.9])|

|Voltages on DDS 9|Total failing rate|
|-----------------|------------------|
|3.075, 1.608     |0.00%             |
|3.107, 1.641     |1.21%             |
|3.195, 1.670     |43.51%            |

## 2014/12/23

### Test condition

`0x80` write with each step `0x10000` repeated `0x100` times on each channel.

### Before any change to the hardware

Connection one thin wires for both breakout:
```
Total failing rate: 0.00%
```

Two thin wires for breakout 2, one thick wire for breakout 1:
```
Failing rate for DDS 3: 0.01+-0.01%
Failing rate for DDS 4: 0.00+-0.00%
Failing rate for DDS 5: 0.02+-0.01%
Failing rate for DDS 7: 6.66+-0.30%
Failing rate for DDS 8: 1.00+-0.07%
Total failing rate: 0.37%
```

Two thin wires for breakout 2, one thick wire for breakout 1 (10% trim up):
```
Failing rate for DDS 0: 99.22+-0.00%
Failing rate for DDS 1: 99.22+-0.00%
Failing rate for DDS 2: 99.22+-0.00%
Failing rate for DDS 3: 99.22+-0.00%
Failing rate for DDS 4: 99.22+-0.00%
Failing rate for DDS 5: 99.22+-0.00%
Failing rate for DDS 6: 99.22+-0.00%
Failing rate for DDS 7: 99.22+-0.00%
Failing rate for DDS 8: 99.22+-0.00%
Failing rate for DDS 9: 99.22+-0.00%
Failing rate for DDS 11: 60.87+-0.32%
Failing rate for DDS 12: 59.91+-0.32%
Failing rate for DDS 13: 99.69+-0.02%
Failing rate for DDS 14: 99.23+-0.01%
Failing rate for DDS 15: 99.22+-0.00%
Failing rate for DDS 16: 99.22+-0.00%
Failing rate for DDS 17: 99.74+-0.02%
Failing rate for DDS 18: 99.74+-0.02%
Failing rate for DDS 19: 99.65+-0.02%
Failing rate for DDS 20: 99.57+-0.02%
Failing rate for DDS 21: 87.57+-0.40%
Total failing rate: 95.08%
```
