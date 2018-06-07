EESchema Schematic File Version 2
LIBS:power
LIBS:device
LIBS:switches
LIBS:relays
LIBS:motors
LIBS:transistors
LIBS:conn
LIBS:linear
LIBS:regul
LIBS:74xx
LIBS:cmos4000
LIBS:adc-dac
LIBS:memory
LIBS:xilinx
LIBS:microcontrollers
LIBS:dsp
LIBS:microchip
LIBS:analog_switches
LIBS:motorola
LIBS:texas
LIBS:intel
LIBS:audio
LIBS:interface
LIBS:digital-audio
LIBS:philips
LIBS:display
LIBS:cypress
LIBS:siliconi
LIBS:opto
LIBS:atmel
LIBS:contrib
LIBS:valves
LIBS:regulator
LIBS:regulator-cache
EELAYER 25 0
EELAYER END
$Descr A4 11693 8268
encoding utf-8
Sheet 1 1
Title ""
Date ""
Rev ""
Comp ""
Comment1 ""
Comment2 ""
Comment3 ""
Comment4 ""
$EndDescr
$Comp
L LT3089-TSSOP U1
U 1 1 5B1970B7
P 5250 3050
F 0 "U1" H 5000 3550 50  0000 C CNN
F 1 "LT3089-TSSOP" H 5100 3550 50  0000 L CNN
F 2 "Housings_SSOP:TSSOP-16-1EP_4.4x5mm_Pitch0.65mm" H 5200 3650 50  0001 C CNN
F 3 "" H 5200 3650 50  0001 C CNN
	1    5250 3050
	-1   0    0    1   
$EndComp
$Comp
L LT3090-MSOP U2
U 1 1 5B197128
P 5250 4500
F 0 "U2" H 5050 4900 50  0000 C CNN
F 1 "LT3090-MSOP" H 5250 4900 50  0000 L CNN
F 2 "Housings_SSOP:MSOP-12-1EP_3x4mm_Pitch0.65mm" H 5200 5000 50  0001 C CNN
F 3 "" H 5200 5100 50  0001 C CNN
	1    5250 4500
	1    0    0    -1  
$EndComp
$Comp
L Conn_01x03 J1
U 1 1 5B1971ED
P 1950 3850
F 0 "J1" H 1950 4050 50  0000 C CNN
F 1 "Conn_01x03" H 1950 3650 50  0000 C CNN
F 2 "Connectors2:3-pins_150mil" H 1950 3850 50  0001 C CNN
F 3 "" H 1950 3850 50  0001 C CNN
	1    1950 3850
	-1   0    0    -1  
$EndComp
$Comp
L Conn_01x03 J2
U 1 1 5B1972F9
P 9800 3900
F 0 "J2" H 9800 4100 50  0000 C CNN
F 1 "Conn_01x03" H 9800 3700 50  0000 C CNN
F 2 "Connectors2:3-pins_150mil" H 9800 3900 50  0001 C CNN
F 3 "" H 9800 3900 50  0001 C CNN
	1    9800 3900
	1    0    0    -1  
$EndComp
Text GLabel 2250 3750 2    60   Input ~ 0
VIN+
Text GLabel 2250 3950 2    60   Input ~ 0
VIN-
$Comp
L GND #PWR01
U 1 1 5B197434
P 2650 3900
F 0 "#PWR01" H 2650 3650 50  0001 C CNN
F 1 "GND" H 2650 3750 50  0000 C CNN
F 2 "" H 2650 3900 50  0001 C CNN
F 3 "" H 2650 3900 50  0001 C CNN
	1    2650 3900
	1    0    0    -1  
$EndComp
Text Label 5700 2700 0    60   ~ 0
RVOUT+
Text Label 5250 2450 0    60   ~ 0
RVOUT+
Text Label 4800 2800 2    60   ~ 0
RVOUT+
Text Label 4800 3400 3    60   ~ 0
RVOUT+
Text Label 5750 3200 3    60   ~ 0
RVOUT+
NoConn ~ 5600 2900
$Comp
L C C10
U 1 1 5B197A4B
P 4600 3250
F 0 "C10" H 4625 3350 50  0000 L CNN
F 1 "0.1u" H 4625 3150 50  0000 L CNN
F 2 "Capacitors_SMD:C_0805" H 4638 3100 50  0001 C CNN
F 3 "" H 4600 3250 50  0001 C CNN
	1    4600 3250
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR02
U 1 1 5B197AA9
P 4600 3500
F 0 "#PWR02" H 4600 3250 50  0001 C CNN
F 1 "GND" H 4600 3350 50  0000 C CNN
F 2 "" H 4600 3500 50  0001 C CNN
F 3 "" H 4600 3500 50  0001 C CNN
	1    4600 3500
	1    0    0    -1  
$EndComp
$Comp
L R R2
U 1 1 5B197BF2
P 3950 3000
F 0 "R2" V 4030 3000 50  0000 C CNN
F 1 "0.1" V 3950 3000 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 3880 3000 50  0001 C CNN
F 3 "" H 3950 3000 50  0001 C CNN
	1    3950 3000
	0    1    1    0   
$EndComp
$Comp
L C C8
U 1 1 5B197C5F
P 4400 3250
F 0 "C8" H 4425 3350 50  0000 L CNN
F 1 "1.5u" H 4425 3150 50  0000 L CNN
F 2 "Capacitors_SMD:C_0805" H 4438 3100 50  0001 C CNN
F 3 "" H 4400 3250 50  0001 C CNN
	1    4400 3250
	1    0    0    -1  
$EndComp
$Comp
L C C6
U 1 1 5B197C88
P 4200 3250
F 0 "C6" H 4225 3350 50  0000 L CNN
F 1 "10u" H 4225 3150 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 4238 3100 50  0001 C CNN
F 3 "" H 4200 3250 50  0001 C CNN
	1    4200 3250
	1    0    0    -1  
$EndComp
$Comp
L CP C2
U 1 1 5B197EE8
P 3450 3250
F 0 "C2" H 3475 3350 50  0000 L CNN
F 1 "100u" H 3250 3150 50  0000 L CNN
F 2 "Capacitors_THT:CP_Radial_D10.0mm_P5.00mm" H 3488 3100 50  0001 C CNN
F 3 "" H 3450 3250 50  0001 C CNN
	1    3450 3250
	1    0    0    -1  
$EndComp
$Comp
L C C4
U 1 1 5B197F2C
P 3700 3250
F 0 "C4" H 3725 3350 50  0000 L CNN
F 1 "10u" H 3725 3150 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 3738 3100 50  0001 C CNN
F 3 "" H 3700 3250 50  0001 C CNN
	1    3700 3250
	1    0    0    -1  
$EndComp
Text GLabel 3300 3000 0    60   Input ~ 0
VIN+
$Comp
L R R3
U 1 1 5B19867A
P 6000 3050
F 0 "R3" V 6080 3050 50  0000 C CNN
F 1 "232k" V 6000 3050 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 5930 3050 50  0001 C CNN
F 3 "" H 6000 3050 50  0001 C CNN
	1    6000 3050
	-1   0    0    1   
$EndComp
$Comp
L POT RV1
U 1 1 5B1986E9
P 6000 3450
F 0 "RV1" V 5825 3450 50  0000 C CNN
F 1 "20k" V 5900 3450 50  0000 C CNN
F 2 "Potentiometers:Potentiometer_Trimmer_Bourns_3269W" H 6000 3450 50  0001 C CNN
F 3 "" H 6000 3450 50  0001 C CNN
	1    6000 3450
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR03
U 1 1 5B19887F
P 6000 3700
F 0 "#PWR03" H 6000 3450 50  0001 C CNN
F 1 "GND" H 6000 3550 50  0000 C CNN
F 2 "" H 6000 3700 50  0001 C CNN
F 3 "" H 6000 3700 50  0001 C CNN
	1    6000 3700
	1    0    0    -1  
$EndComp
$Comp
L C C11
U 1 1 5B19897B
P 6300 3450
F 0 "C11" H 6325 3550 50  0000 L CNN
F 1 "10u" H 6325 3350 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 6338 3300 50  0001 C CNN
F 3 "" H 6300 3450 50  0001 C CNN
	1    6300 3450
	1    0    0    -1  
$EndComp
$Comp
L C C13
U 1 1 5B198E04
P 6500 3250
F 0 "C13" H 6525 3350 50  0000 L CNN
F 1 "0.1u" H 6525 3150 50  0000 L CNN
F 2 "Capacitors_SMD:C_0805" H 6538 3100 50  0001 C CNN
F 3 "" H 6500 3250 50  0001 C CNN
	1    6500 3250
	1    0    0    -1  
$EndComp
$Comp
L C C15
U 1 1 5B198E4E
P 6700 3050
F 0 "C15" H 6725 3150 50  0000 L CNN
F 1 "1.5u" H 6725 2950 50  0000 L CNN
F 2 "Capacitors_SMD:C_0805" H 6738 2900 50  0001 C CNN
F 3 "" H 6700 3050 50  0001 C CNN
	1    6700 3050
	1    0    0    -1  
$EndComp
Text Label 7250 2700 2    60   ~ 0
RVOUT+
$Comp
L C C18
U 1 1 5B19958E
P 7600 2900
F 0 "C18" H 7625 3000 50  0000 L CNN
F 1 "10u" H 7625 2800 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 7638 2750 50  0001 C CNN
F 3 "" H 7600 2900 50  0001 C CNN
	1    7600 2900
	1    0    0    -1  
$EndComp
$Comp
L C C20
U 1 1 5B199644
P 7800 2900
F 0 "C20" H 7825 3000 50  0000 L CNN
F 1 "1.5u" H 7825 2800 50  0000 L CNN
F 2 "Capacitors_SMD:C_0805" H 7838 2750 50  0001 C CNN
F 3 "" H 7800 2900 50  0001 C CNN
	1    7800 2900
	1    0    0    -1  
$EndComp
$Comp
L C C22
U 1 1 5B1996D0
P 8000 2900
F 0 "C22" H 8025 3000 50  0000 L CNN
F 1 "0.1u" H 8025 2800 50  0000 L CNN
F 2 "Capacitors_SMD:C_0805" H 8038 2750 50  0001 C CNN
F 3 "" H 8000 2900 50  0001 C CNN
	1    8000 2900
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR04
U 1 1 5B199830
P 7600 3150
F 0 "#PWR04" H 7600 2900 50  0001 C CNN
F 1 "GND" H 7600 3000 50  0000 C CNN
F 2 "" H 7600 3150 50  0001 C CNN
F 3 "" H 7600 3150 50  0001 C CNN
	1    7600 3150
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR05
U 1 1 5B19A16D
P 4800 4750
F 0 "#PWR05" H 4800 4500 50  0001 C CNN
F 1 "GND" H 4800 4600 50  0000 C CNN
F 2 "" H 4800 4750 50  0001 C CNN
F 3 "" H 4800 4750 50  0001 C CNN
	1    4800 4750
	0    1    1    0   
$EndComp
$Comp
L GND #PWR06
U 1 1 5B19A44B
P 5700 4550
F 0 "#PWR06" H 5700 4300 50  0001 C CNN
F 1 "GND" H 5700 4400 50  0000 C CNN
F 2 "" H 5700 4550 50  0001 C CNN
F 3 "" H 5700 4550 50  0001 C CNN
	1    5700 4550
	0    -1   -1   0   
$EndComp
Text Label 4800 4150 2    60   ~ 0
RVIN-
Text Label 5750 4750 3    60   ~ 0
RVIN-
Text Label 5250 5050 0    60   ~ 0
RVIN-
$Comp
L C C9
U 1 1 5B19B168
P 4400 4700
F 0 "C9" H 4425 4800 50  0000 L CNN
F 1 "0.1u" H 4425 4600 50  0000 L CNN
F 2 "Capacitors_SMD:C_0805" H 4438 4550 50  0001 C CNN
F 3 "" H 4400 4700 50  0001 C CNN
	1    4400 4700
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR07
U 1 1 5B19B16E
P 4400 4950
F 0 "#PWR07" H 4400 4700 50  0001 C CNN
F 1 "GND" H 4400 4800 50  0000 C CNN
F 2 "" H 4400 4950 50  0001 C CNN
F 3 "" H 4400 4950 50  0001 C CNN
	1    4400 4950
	1    0    0    -1  
$EndComp
$Comp
L R R1
U 1 1 5B19B174
P 3750 4450
F 0 "R1" V 3830 4450 50  0000 C CNN
F 1 "0.1" V 3750 4450 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 3680 4450 50  0001 C CNN
F 3 "" H 3750 4450 50  0001 C CNN
	1    3750 4450
	0    1    1    0   
$EndComp
$Comp
L C C7
U 1 1 5B19B17A
P 4200 4700
F 0 "C7" H 4225 4800 50  0000 L CNN
F 1 "1.5u" H 4225 4600 50  0000 L CNN
F 2 "Capacitors_SMD:C_0805" H 4238 4550 50  0001 C CNN
F 3 "" H 4200 4700 50  0001 C CNN
	1    4200 4700
	1    0    0    -1  
$EndComp
$Comp
L C C5
U 1 1 5B19B180
P 4000 4700
F 0 "C5" H 4025 4800 50  0000 L CNN
F 1 "10u" H 4025 4600 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 4038 4550 50  0001 C CNN
F 3 "" H 4000 4700 50  0001 C CNN
	1    4000 4700
	1    0    0    -1  
$EndComp
$Comp
L CP C1
U 1 1 5B19B186
P 3250 4700
F 0 "C1" H 3275 4800 50  0000 L CNN
F 1 "100u" H 3050 4600 50  0000 L CNN
F 2 "Capacitors_THT:CP_Radial_D10.0mm_P5.00mm" H 3288 4550 50  0001 C CNN
F 3 "" H 3250 4700 50  0001 C CNN
	1    3250 4700
	1    0    0    1   
$EndComp
$Comp
L C C3
U 1 1 5B19B18C
P 3500 4700
F 0 "C3" H 3525 4800 50  0000 L CNN
F 1 "10u" H 3525 4600 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 3538 4550 50  0001 C CNN
F 3 "" H 3500 4700 50  0001 C CNN
	1    3500 4700
	1    0    0    -1  
$EndComp
Text GLabel 3100 4450 0    60   Input ~ 0
VIN-
$Comp
L R R4
U 1 1 5B19B927
P 6050 4900
F 0 "R4" V 6130 4900 50  0000 C CNN
F 1 "232k" V 6050 4900 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 5980 4900 50  0001 C CNN
F 3 "" H 6050 4900 50  0001 C CNN
	1    6050 4900
	-1   0    0    1   
$EndComp
$Comp
L POT RV2
U 1 1 5B19B92D
P 6050 5300
F 0 "RV2" V 5875 5300 50  0000 C CNN
F 1 "20k" V 5950 5300 50  0000 C CNN
F 2 "Potentiometers:Potentiometer_Trimmer_Bourns_3269W" H 6050 5300 50  0001 C CNN
F 3 "" H 6050 5300 50  0001 C CNN
	1    6050 5300
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR08
U 1 1 5B19B933
P 6050 5550
F 0 "#PWR08" H 6050 5300 50  0001 C CNN
F 1 "GND" H 6050 5400 50  0000 C CNN
F 2 "" H 6050 5550 50  0001 C CNN
F 3 "" H 6050 5550 50  0001 C CNN
	1    6050 5550
	1    0    0    -1  
$EndComp
$Comp
L C C12
U 1 1 5B19B939
P 6350 5300
F 0 "C12" H 6375 5400 50  0000 L CNN
F 1 "10u" H 6375 5200 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 6388 5150 50  0001 C CNN
F 3 "" H 6350 5300 50  0001 C CNN
	1    6350 5300
	1    0    0    -1  
$EndComp
$Comp
L C C14
U 1 1 5B19B93F
P 6550 5100
F 0 "C14" H 6575 5200 50  0000 L CNN
F 1 "0.1u" H 6575 5000 50  0000 L CNN
F 2 "Capacitors_SMD:C_0805" H 6588 4950 50  0001 C CNN
F 3 "" H 6550 5100 50  0001 C CNN
	1    6550 5100
	1    0    0    -1  
$EndComp
$Comp
L C C16
U 1 1 5B19B945
P 6750 4900
F 0 "C16" H 6775 5000 50  0000 L CNN
F 1 "1.5u" H 6775 4800 50  0000 L CNN
F 2 "Capacitors_SMD:C_0805" H 6788 4750 50  0001 C CNN
F 3 "" H 6750 4900 50  0001 C CNN
	1    6750 4900
	1    0    0    -1  
$EndComp
Text Label 5800 4350 0    60   ~ 0
RVOUT-
$Comp
L C C17
U 1 1 5B19BEF0
P 7550 4500
F 0 "C17" H 7575 4600 50  0000 L CNN
F 1 "10u" H 7575 4400 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 7588 4350 50  0001 C CNN
F 3 "" H 7550 4500 50  0001 C CNN
	1    7550 4500
	1    0    0    -1  
$EndComp
$Comp
L C C19
U 1 1 5B19BEF6
P 7750 4500
F 0 "C19" H 7775 4600 50  0000 L CNN
F 1 "1.5u" H 7775 4400 50  0000 L CNN
F 2 "Capacitors_SMD:C_0805" H 7788 4350 50  0001 C CNN
F 3 "" H 7750 4500 50  0001 C CNN
	1    7750 4500
	1    0    0    -1  
$EndComp
$Comp
L C C21
U 1 1 5B19BEFC
P 7950 4500
F 0 "C21" H 7975 4600 50  0000 L CNN
F 1 "0.1u" H 7975 4400 50  0000 L CNN
F 2 "Capacitors_SMD:C_0805" H 7988 4350 50  0001 C CNN
F 3 "" H 7950 4500 50  0001 C CNN
	1    7950 4500
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR09
U 1 1 5B19BF08
P 7550 4750
F 0 "#PWR09" H 7550 4500 50  0001 C CNN
F 1 "GND" H 7550 4600 50  0000 C CNN
F 2 "" H 7550 4750 50  0001 C CNN
F 3 "" H 7550 4750 50  0001 C CNN
	1    7550 4750
	1    0    0    -1  
$EndComp
Text GLabel 9500 3800 0    60   Input ~ 0
VOUT+
Text GLabel 8100 2700 2    60   Input ~ 0
VOUT+
Text GLabel 9500 4000 0    60   Input ~ 0
VOUT-
Text GLabel 8050 4300 2    60   Input ~ 0
VOUT-
$Comp
L GND #PWR010
U 1 1 5B19C8A9
P 9050 3950
F 0 "#PWR010" H 9050 3700 50  0001 C CNN
F 1 "GND" H 9050 3800 50  0000 C CNN
F 2 "" H 9050 3950 50  0001 C CNN
F 3 "" H 9050 3950 50  0001 C CNN
	1    9050 3950
	1    0    0    -1  
$EndComp
Wire Wire Line
	2150 3750 2250 3750
Wire Wire Line
	2150 3950 2250 3950
Wire Wire Line
	2150 3850 2650 3850
Wire Wire Line
	2650 3850 2650 3900
Wire Wire Line
	5600 2700 5700 2700
Wire Wire Line
	5250 2450 5250 2500
Wire Wire Line
	4800 2800 4900 2800
Wire Wire Line
	4850 2700 4850 2900
Wire Wire Line
	4850 2700 4900 2700
Connection ~ 4850 2800
Wire Wire Line
	4850 2900 4900 2900
Wire Wire Line
	4800 3400 4900 3400
Wire Wire Line
	5600 3200 5750 3200
Wire Wire Line
	5600 3000 5650 3000
Wire Wire Line
	5650 3000 5650 3400
Connection ~ 5650 3200
Wire Wire Line
	5600 3100 5650 3100
Connection ~ 5650 3100
Wire Wire Line
	5650 3300 5600 3300
Wire Wire Line
	5650 3400 5600 3400
Connection ~ 5650 3300
Wire Wire Line
	4100 3000 4900 3000
Wire Wire Line
	4850 3000 4850 3300
Wire Wire Line
	4850 3300 4900 3300
Wire Wire Line
	4900 3200 4850 3200
Connection ~ 4850 3200
Wire Wire Line
	4900 3100 4850 3100
Connection ~ 4850 3100
Wire Wire Line
	4600 3500 4600 3400
Wire Wire Line
	4600 3100 4600 3000
Connection ~ 4850 3000
Wire Wire Line
	4200 3400 4200 3450
Wire Wire Line
	3450 3450 4600 3450
Connection ~ 4600 3450
Wire Wire Line
	4400 3400 4400 3450
Connection ~ 4400 3450
Connection ~ 4600 3000
Wire Wire Line
	4400 3100 4400 3000
Connection ~ 4400 3000
Wire Wire Line
	4200 3100 4200 3000
Connection ~ 4200 3000
Wire Wire Line
	3300 3000 3800 3000
Wire Wire Line
	3700 3000 3700 3100
Wire Wire Line
	3450 3000 3450 3100
Connection ~ 3700 3000
Wire Wire Line
	3450 3400 3450 3450
Connection ~ 4200 3450
Wire Wire Line
	3700 3400 3700 3450
Connection ~ 3700 3450
Connection ~ 3450 3000
Wire Wire Line
	5600 2800 6700 2800
Wire Wire Line
	6000 2800 6000 2900
Wire Wire Line
	6000 3200 6000 3300
Connection ~ 6000 3250
Wire Wire Line
	6000 3600 6000 3700
Wire Wire Line
	6150 3650 6150 3450
Connection ~ 6000 3650
Connection ~ 6150 3650
Wire Wire Line
	6000 3250 6300 3250
Wire Wire Line
	6300 3250 6300 3300
Wire Wire Line
	6300 3650 6300 3600
Wire Wire Line
	6000 3650 6700 3650
Wire Wire Line
	6500 2800 6500 3100
Connection ~ 6000 2800
Wire Wire Line
	6500 3650 6500 3400
Connection ~ 6300 3650
Wire Wire Line
	6700 2800 6700 2900
Connection ~ 6500 2800
Wire Wire Line
	6700 3650 6700 3200
Connection ~ 6500 3650
Wire Wire Line
	7250 2700 8100 2700
Wire Wire Line
	7600 2700 7600 2750
Wire Wire Line
	7800 2700 7800 2750
Connection ~ 7600 2700
Wire Wire Line
	8000 2700 8000 2750
Connection ~ 7800 2700
Wire Wire Line
	7600 3050 7600 3150
Wire Wire Line
	7800 3100 7800 3050
Connection ~ 7600 3100
Wire Wire Line
	8000 3100 8000 3050
Connection ~ 7800 3100
Connection ~ 8000 2700
Wire Wire Line
	4850 4250 4800 4250
Wire Wire Line
	4800 4150 4800 4650
Wire Wire Line
	4800 4350 4850 4350
Wire Wire Line
	3900 4450 4850 4450
Connection ~ 4800 4350
Wire Wire Line
	4800 4550 4850 4550
Connection ~ 4800 4450
Wire Wire Line
	4800 4650 4850 4650
Connection ~ 4800 4550
Wire Wire Line
	4800 4750 4850 4750
Wire Wire Line
	5700 4550 5650 4550
Connection ~ 4800 4250
Wire Wire Line
	5650 4750 5750 4750
Wire Wire Line
	5700 4450 5650 4450
Wire Wire Line
	5700 4250 5700 4450
Wire Wire Line
	5650 4350 5800 4350
Wire Wire Line
	5700 4250 5650 4250
Connection ~ 5700 4350
Wire Wire Line
	4400 4950 4400 4850
Wire Wire Line
	4400 4550 4400 4450
Wire Wire Line
	4000 4850 4000 4900
Wire Wire Line
	3250 4900 4400 4900
Connection ~ 4400 4900
Wire Wire Line
	4200 4850 4200 4900
Connection ~ 4200 4900
Connection ~ 4400 4450
Wire Wire Line
	4200 4550 4200 4450
Connection ~ 4200 4450
Wire Wire Line
	4000 4550 4000 4450
Connection ~ 4000 4450
Wire Wire Line
	3100 4450 3600 4450
Wire Wire Line
	3500 4450 3500 4550
Wire Wire Line
	3250 4450 3250 4550
Connection ~ 3500 4450
Wire Wire Line
	3250 4850 3250 4900
Connection ~ 4000 4900
Wire Wire Line
	3500 4850 3500 4900
Connection ~ 3500 4900
Connection ~ 3250 4450
Wire Wire Line
	5650 4650 6750 4650
Wire Wire Line
	6050 4650 6050 4750
Wire Wire Line
	6050 5050 6050 5150
Connection ~ 6050 5100
Wire Wire Line
	6050 5450 6050 5550
Wire Wire Line
	6200 5500 6200 5300
Connection ~ 6050 5500
Connection ~ 6200 5500
Wire Wire Line
	6050 5100 6350 5100
Wire Wire Line
	6350 5100 6350 5150
Wire Wire Line
	6350 5500 6350 5450
Wire Wire Line
	6050 5500 6750 5500
Wire Wire Line
	6550 4650 6550 4950
Connection ~ 6050 4650
Wire Wire Line
	6550 5500 6550 5250
Connection ~ 6350 5500
Wire Wire Line
	6750 4650 6750 4750
Connection ~ 6550 4650
Wire Wire Line
	6750 5500 6750 5050
Connection ~ 6550 5500
Wire Wire Line
	7250 4300 8050 4300
Wire Wire Line
	7550 4300 7550 4350
Wire Wire Line
	7750 4300 7750 4350
Connection ~ 7550 4300
Wire Wire Line
	7950 4300 7950 4350
Connection ~ 7750 4300
Wire Wire Line
	7550 4650 7550 4750
Wire Wire Line
	7350 4700 7950 4700
Wire Wire Line
	7750 4700 7750 4650
Connection ~ 7550 4700
Wire Wire Line
	7950 4700 7950 4650
Connection ~ 7750 4700
Connection ~ 7950 4300
Wire Wire Line
	9050 3950 9050 3900
Wire Wire Line
	9050 3900 9600 3900
Wire Wire Line
	9600 3800 9500 3800
Wire Wire Line
	9500 4000 9600 4000
$Comp
L R R5
U 1 1 5B19D2D6
P 7350 2900
F 0 "R5" V 7430 2900 50  0000 C CNN
F 1 "3.3k" V 7350 2900 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 7280 2900 50  0001 C CNN
F 3 "" H 7350 2900 50  0001 C CNN
	1    7350 2900
	-1   0    0    1   
$EndComp
Wire Wire Line
	7350 2700 7350 2750
Connection ~ 7350 2700
Wire Wire Line
	7350 3050 7350 3100
Text Label 7250 4300 2    60   ~ 0
RVOUT-
$Comp
L R R6
U 1 1 5B19DBB7
P 7350 4500
F 0 "R6" V 7430 4500 50  0000 C CNN
F 1 "3.3k" V 7350 4500 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 7280 4500 50  0001 C CNN
F 3 "" H 7350 4500 50  0001 C CNN
	1    7350 4500
	-1   0    0    1   
$EndComp
Wire Wire Line
	7350 4350 7350 4300
Connection ~ 7350 4300
Wire Wire Line
	7350 4650 7350 4700
Wire Wire Line
	7350 3100 8000 3100
Wire Wire Line
	5250 5050 5250 4950
$EndSCHEMATC
