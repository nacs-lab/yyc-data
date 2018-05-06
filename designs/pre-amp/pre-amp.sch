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
LIBS:pre-amp-cache
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
NoConn ~ 2000 1400
NoConn ~ 2000 1600
$Comp
L Conn_Coaxial J1
U 1 1 5AEB2F69
P 750 1300
F 0 "J1" H 760 1420 50  0000 C CNN
F 1 "Conn_Coaxial" V 865 1300 50  0000 C CNN
F 2 "Connectors_TE-Connectivity:BNC_Socket_TYCO-AMP_LargePads" H 750 1300 50  0001 C CNN
F 3 "" H 750 1300 50  0001 C CNN
	1    750  1300
	0    -1   -1   0   
$EndComp
$Comp
L GND #PWR13
U 1 1 5AEB2FF5
P 2500 1900
F 0 "#PWR13" H 2500 1650 50  0001 C CNN
F 1 "GND" H 2500 1750 50  0000 C CNN
F 2 "" H 2500 1900 50  0001 C CNN
F 3 "" H 2500 1900 50  0001 C CNN
	1    2500 1900
	1    0    0    -1  
$EndComp
$Comp
L C C1
U 1 1 5AEB3124
P 1350 2700
F 0 "C1" H 1375 2800 50  0000 L CNN
F 1 "0.1u" H 1375 2600 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 1388 2550 50  0001 C CNN
F 3 "" H 1350 2700 50  0001 C CNN
	1    1350 2700
	1    0    0    -1  
$EndComp
$Comp
L C C3
U 1 1 5AEB3174
P 1600 2700
F 0 "C3" H 1625 2800 50  0000 L CNN
F 1 "10u" H 1625 2600 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 1638 2550 50  0001 C CNN
F 3 "" H 1600 2700 50  0001 C CNN
	1    1600 2700
	1    0    0    -1  
$EndComp
$Comp
L -15V #PWR8
U 1 1 5AEB31DF
P 2300 1900
F 0 "#PWR8" H 2300 2000 50  0001 C CNN
F 1 "-15V" H 2300 2050 50  0000 C CNN
F 2 "" H 2300 1900 50  0001 C CNN
F 3 "" H 2300 1900 50  0001 C CNN
	1    2300 1900
	-1   0    0    1   
$EndComp
$Comp
L GND #PWR3
U 1 1 5AEB3300
P 1350 2950
F 0 "#PWR3" H 1350 2700 50  0001 C CNN
F 1 "GND" H 1350 2800 50  0000 C CNN
F 2 "" H 1350 2950 50  0001 C CNN
F 3 "" H 1350 2950 50  0001 C CNN
	1    1350 2950
	1    0    0    -1  
$EndComp
$Comp
L -15V #PWR2
U 1 1 5AEB33C1
P 1350 2450
F 0 "#PWR2" H 1350 2550 50  0001 C CNN
F 1 "-15V" H 1350 2600 50  0000 C CNN
F 2 "" H 1350 2450 50  0001 C CNN
F 3 "" H 1350 2450 50  0001 C CNN
	1    1350 2450
	1    0    0    -1  
$EndComp
$Comp
L +15V #PWR7
U 1 1 5AEB347E
P 2300 1100
F 0 "#PWR7" H 2300 950 50  0001 C CNN
F 1 "+15V" H 2300 1240 50  0000 C CNN
F 2 "" H 2300 1100 50  0001 C CNN
F 3 "" H 2300 1100 50  0001 C CNN
	1    2300 1100
	1    0    0    -1  
$EndComp
$Comp
L +15V #PWR6
U 1 1 5AEB34BD
P 1850 2450
F 0 "#PWR6" H 1850 2300 50  0001 C CNN
F 1 "+15V" H 1850 2590 50  0000 C CNN
F 2 "" H 1850 2450 50  0001 C CNN
F 3 "" H 1850 2450 50  0001 C CNN
	1    1850 2450
	1    0    0    -1  
$EndComp
$Comp
L C C4
U 1 1 5AEB3510
P 1850 2700
F 0 "C4" H 1875 2800 50  0000 L CNN
F 1 "0.1u" H 1875 2600 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 1888 2550 50  0001 C CNN
F 3 "" H 1850 2700 50  0001 C CNN
	1    1850 2700
	1    0    0    -1  
$EndComp
$Comp
L C C5
U 1 1 5AEB3546
P 2100 2700
F 0 "C5" H 2125 2800 50  0000 L CNN
F 1 "10u" H 2125 2600 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 2138 2550 50  0001 C CNN
F 3 "" H 2100 2700 50  0001 C CNN
	1    2100 2700
	1    0    0    -1  
$EndComp
Wire Wire Line
	2300 1800 2300 1900
Wire Wire Line
	2500 1800 2500 1900
Wire Wire Line
	1350 2850 1350 2950
Wire Wire Line
	1600 2850 1600 2900
Wire Wire Line
	1350 2900 2100 2900
Connection ~ 1350 2900
Wire Wire Line
	1600 2550 1600 2500
Wire Wire Line
	1600 2500 1350 2500
Wire Wire Line
	1350 2450 1350 2550
Connection ~ 1350 2500
Wire Wire Line
	2300 1100 2300 1200
Wire Wire Line
	1850 2450 1850 2550
Wire Wire Line
	1850 2900 1850 2850
Connection ~ 1600 2900
Wire Wire Line
	1850 2500 2100 2500
Wire Wire Line
	2100 2500 2100 2550
Connection ~ 1850 2500
Wire Wire Line
	2100 2900 2100 2850
Connection ~ 1850 2900
$Comp
L R R1
U 1 1 5AEB39CD
P 1200 1300
F 0 "R1" V 1280 1300 50  0000 C CNN
F 1 "1k" V 1200 1300 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 1130 1300 50  0001 C CNN
F 3 "" H 1200 1300 50  0001 C CNN
	1    1200 1300
	0    1    1    0   
$EndComp
$Comp
L C C2
U 1 1 5AEB3A3A
P 1450 1500
F 0 "C2" H 1475 1600 50  0000 L CNN
F 1 "0.02u" H 1475 1400 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 1488 1350 50  0001 C CNN
F 3 "" H 1450 1500 50  0001 C CNN
	1    1450 1500
	1    0    0    -1  
$EndComp
$Comp
L R R2
U 1 1 5AEB3AAC
P 1200 1700
F 0 "R2" V 1280 1700 50  0000 C CNN
F 1 "1k" V 1200 1700 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 1130 1700 50  0001 C CNN
F 3 "" H 1200 1700 50  0001 C CNN
	1    1200 1700
	0    1    1    0   
$EndComp
Wire Wire Line
	750  1450 750  1700
Wire Wire Line
	750  1700 1050 1700
Wire Wire Line
	950  1300 1050 1300
$Comp
L R R4
U 1 1 5AEB3C27
P 1850 1500
F 0 "R4" V 1930 1500 50  0000 C CNN
F 1 "20k" V 1850 1500 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 1780 1500 50  0001 C CNN
F 3 "" H 1850 1500 50  0001 C CNN
	1    1850 1500
	-1   0    0    1   
$EndComp
Wire Wire Line
	1350 1300 2000 1300
Wire Wire Line
	2000 1700 1350 1700
Wire Wire Line
	1450 1650 1450 1700
Connection ~ 1450 1700
Wire Wire Line
	1450 1350 1450 1300
Connection ~ 1450 1300
Wire Wire Line
	1850 1650 1850 1700
Connection ~ 1850 1700
Wire Wire Line
	1850 1350 1850 1300
Connection ~ 1850 1300
$Comp
L R R5
U 1 1 5AEB3FD7
P 3100 1500
F 0 "R5" V 3180 1500 50  0000 C CNN
F 1 "R" V 3100 1500 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 3030 1500 50  0001 C CNN
F 3 "" H 3100 1500 50  0001 C CNN
	1    3100 1500
	0    1    1    0   
$EndComp
$Comp
L C C12
U 1 1 5AEB405F
P 3550 1750
F 0 "C12" H 3575 1850 50  0000 L CNN
F 1 "C" H 3575 1650 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 3588 1600 50  0001 C CNN
F 3 "" H 3550 1750 50  0001 C CNN
	1    3550 1750
	1    0    0    -1  
$EndComp
$Comp
L R R7
U 1 1 5AEB40A1
P 3800 1500
F 0 "R7" V 3880 1500 50  0000 C CNN
F 1 "R" V 3800 1500 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 3730 1500 50  0001 C CNN
F 3 "" H 3800 1500 50  0001 C CNN
	1    3800 1500
	0    1    1    0   
$EndComp
$Comp
L C C13
U 1 1 5AEB41BC
P 4500 900
F 0 "C13" H 4525 1000 50  0000 L CNN
F 1 "C" H 4525 800 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 4538 750 50  0001 C CNN
F 3 "" H 4500 900 50  0001 C CNN
	1    4500 900 
	0    1    1    0   
$EndComp
$Comp
L R R6
U 1 1 5AEB4241
P 3800 750
F 0 "R6" V 3880 750 50  0000 C CNN
F 1 "R" V 3800 750 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 3730 750 50  0001 C CNN
F 3 "" H 3800 750 50  0001 C CNN
	1    3800 750 
	0    1    1    0   
$EndComp
$Comp
L GND #PWR16
U 1 1 5AEB4363
P 3550 2000
F 0 "#PWR16" H 3550 1750 50  0001 C CNN
F 1 "GND" H 3550 1850 50  0000 C CNN
F 2 "" H 3550 2000 50  0001 C CNN
F 3 "" H 3550 2000 50  0001 C CNN
	1    3550 2000
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR17
U 1 1 5AEB439E
P 4000 2000
F 0 "#PWR17" H 4000 1750 50  0001 C CNN
F 1 "GND" H 4000 1850 50  0000 C CNN
F 2 "" H 4000 2000 50  0001 C CNN
F 3 "" H 4000 2000 50  0001 C CNN
	1    4000 2000
	1    0    0    -1  
$EndComp
Wire Wire Line
	2800 1500 2950 1500
Wire Wire Line
	3250 1500 3650 1500
Wire Wire Line
	3950 1500 4050 1500
Wire Wire Line
	3550 750  3550 1600
Connection ~ 3550 1500
Wire Wire Line
	3550 750  3650 750 
Wire Wire Line
	4650 1600 4850 1600
Wire Wire Line
	4750 750  4750 1600
Wire Wire Line
	4750 900  4650 900 
Wire Wire Line
	4750 750  3950 750 
Connection ~ 4750 900 
Wire Wire Line
	4000 1500 4000 900 
Wire Wire Line
	4000 900  4350 900 
Connection ~ 4000 1500
Wire Wire Line
	4000 2000 4000 1700
Wire Wire Line
	4000 1700 4050 1700
Wire Wire Line
	3550 2000 3550 1900
$Comp
L C C6
U 1 1 5AEB47EF
P 2450 2700
F 0 "C6" H 2475 2800 50  0000 L CNN
F 1 "0.1u" H 2475 2600 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 2488 2550 50  0001 C CNN
F 3 "" H 2450 2700 50  0001 C CNN
	1    2450 2700
	1    0    0    -1  
$EndComp
$Comp
L C C8
U 1 1 5AEB47F5
P 2700 2700
F 0 "C8" H 2725 2800 50  0000 L CNN
F 1 "10u" H 2725 2600 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 2738 2550 50  0001 C CNN
F 3 "" H 2700 2700 50  0001 C CNN
	1    2700 2700
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR10
U 1 1 5AEB47FB
P 2450 2950
F 0 "#PWR10" H 2450 2700 50  0001 C CNN
F 1 "GND" H 2450 2800 50  0000 C CNN
F 2 "" H 2450 2950 50  0001 C CNN
F 3 "" H 2450 2950 50  0001 C CNN
	1    2450 2950
	1    0    0    -1  
$EndComp
$Comp
L -15V #PWR9
U 1 1 5AEB4801
P 2450 2450
F 0 "#PWR9" H 2450 2550 50  0001 C CNN
F 1 "-15V" H 2450 2600 50  0000 C CNN
F 2 "" H 2450 2450 50  0001 C CNN
F 3 "" H 2450 2450 50  0001 C CNN
	1    2450 2450
	1    0    0    -1  
$EndComp
$Comp
L +15V #PWR14
U 1 1 5AEB4807
P 2950 2450
F 0 "#PWR14" H 2950 2300 50  0001 C CNN
F 1 "+15V" H 2950 2590 50  0000 C CNN
F 2 "" H 2950 2450 50  0001 C CNN
F 3 "" H 2950 2450 50  0001 C CNN
	1    2950 2450
	1    0    0    -1  
$EndComp
$Comp
L C C9
U 1 1 5AEB480D
P 2950 2700
F 0 "C9" H 2975 2800 50  0000 L CNN
F 1 "0.1u" H 2975 2600 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 2988 2550 50  0001 C CNN
F 3 "" H 2950 2700 50  0001 C CNN
	1    2950 2700
	1    0    0    -1  
$EndComp
$Comp
L C C11
U 1 1 5AEB4813
P 3200 2700
F 0 "C11" H 3225 2800 50  0000 L CNN
F 1 "10u" H 3225 2600 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 3238 2550 50  0001 C CNN
F 3 "" H 3200 2700 50  0001 C CNN
	1    3200 2700
	1    0    0    -1  
$EndComp
Wire Wire Line
	2450 2850 2450 2950
Wire Wire Line
	2700 2850 2700 2900
Wire Wire Line
	2450 2900 3200 2900
Connection ~ 2450 2900
Wire Wire Line
	2700 2550 2700 2500
Wire Wire Line
	2700 2500 2450 2500
Wire Wire Line
	2450 2450 2450 2550
Connection ~ 2450 2500
Wire Wire Line
	2950 2450 2950 2550
Wire Wire Line
	2950 2900 2950 2850
Connection ~ 2700 2900
Wire Wire Line
	2950 2500 3200 2500
Wire Wire Line
	3200 2500 3200 2550
Connection ~ 2950 2500
Wire Wire Line
	3200 2900 3200 2850
Connection ~ 2950 2900
$Comp
L +15V #PWR19
U 1 1 5AEB4851
P 4250 2000
F 0 "#PWR19" H 4250 1850 50  0001 C CNN
F 1 "+15V" H 4250 2140 50  0000 C CNN
F 2 "" H 4250 2000 50  0001 C CNN
F 3 "" H 4250 2000 50  0001 C CNN
	1    4250 2000
	-1   0    0    1   
$EndComp
$Comp
L -15V #PWR18
U 1 1 5AEB4898
P 4250 1200
F 0 "#PWR18" H 4250 1300 50  0001 C CNN
F 1 "-15V" H 4250 1350 50  0000 C CNN
F 2 "" H 4250 1200 50  0001 C CNN
F 3 "" H 4250 1200 50  0001 C CNN
	1    4250 1200
	1    0    0    -1  
$EndComp
Wire Wire Line
	4250 1200 4250 1300
Wire Wire Line
	4250 1900 4250 2000
$Comp
L R R8
U 1 1 5AEB4F3A
P 5000 1600
F 0 "R8" V 5080 1600 50  0000 C CNN
F 1 "0" V 5000 1600 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 4930 1600 50  0001 C CNN
F 3 "" H 5000 1600 50  0001 C CNN
	1    5000 1600
	0    1    1    0   
$EndComp
$Comp
L C C14
U 1 1 5AEB4FCB
P 5250 1850
F 0 "C14" H 5275 1950 50  0000 L CNN
F 1 "0" H 5275 1750 50  0000 L CNN
F 2 "Capacitors_SMD:C_0603" H 5288 1700 50  0001 C CNN
F 3 "" H 5250 1850 50  0001 C CNN
	1    5250 1850
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR20
U 1 1 5AEB505B
P 5250 2100
F 0 "#PWR20" H 5250 1850 50  0001 C CNN
F 1 "GND" H 5250 1950 50  0000 C CNN
F 2 "" H 5250 2100 50  0001 C CNN
F 3 "" H 5250 2100 50  0001 C CNN
	1    5250 2100
	1    0    0    -1  
$EndComp
Connection ~ 4750 1600
Wire Wire Line
	5150 1600 5450 1600
Wire Wire Line
	5250 1600 5250 1700
Wire Wire Line
	5250 2000 5250 2100
$Comp
L Conn_Coaxial J3
U 1 1 5AEB5216
P 5600 1600
F 0 "J3" H 5610 1720 50  0000 C CNN
F 1 "Conn_Coaxial" V 5715 1600 50  0000 C CNN
F 2 "Connectors_TE-Connectivity:BNC_Socket_TYCO-AMP_LargePads" H 5600 1600 50  0001 C CNN
F 3 "" H 5600 1600 50  0001 C CNN
	1    5600 1600
	1    0    0    -1  
$EndComp
Connection ~ 5250 1600
Wire Wire Line
	5600 1800 5600 2050
Wire Wire Line
	5600 2050 5250 2050
Connection ~ 5250 2050
$Comp
L AD8429 U1
U 1 1 5AEB57AE
P 2400 1500
F 0 "U1" H 2550 1800 50  0000 C CNN
F 1 "AD8429" H 2550 1700 50  0000 C CNN
F 2 "Housings_SOIC:SOIC-8_3.9x4.9mm_Pitch1.27mm" H 2100 1500 50  0001 C CNN
F 3 "" H 2750 1100 50  0001 C CNN
	1    2400 1500
	1    0    0    -1  
$EndComp
$Comp
L Conn_01x03 J2
U 1 1 5AEB6541
P 1750 3550
F 0 "J2" H 1750 3750 50  0000 C CNN
F 1 "Conn_01x03" H 1750 3350 50  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 1750 3550 50  0001 C CNN
F 3 "" H 1750 3550 50  0001 C CNN
	1    1750 3550
	1    0    0    -1  
$EndComp
$Comp
L R R3
U 1 1 5AEB681A
P 1300 3550
F 0 "R3" V 1380 3550 50  0000 C CNN
F 1 "0" V 1300 3550 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 1230 3550 50  0001 C CNN
F 3 "" H 1300 3550 50  0001 C CNN
	1    1300 3550
	0    1    1    0   
$EndComp
$Comp
L GND #PWR1
U 1 1 5AEB68A7
P 1100 3600
F 0 "#PWR1" H 1100 3350 50  0001 C CNN
F 1 "GND" H 1100 3450 50  0000 C CNN
F 2 "" H 1100 3600 50  0001 C CNN
F 3 "" H 1100 3600 50  0001 C CNN
	1    1100 3600
	1    0    0    -1  
$EndComp
$Comp
L -15V #PWR5
U 1 1 5AEB68FD
P 1450 3800
F 0 "#PWR5" H 1450 3900 50  0001 C CNN
F 1 "-15V" H 1450 3950 50  0000 C CNN
F 2 "" H 1450 3800 50  0001 C CNN
F 3 "" H 1450 3800 50  0001 C CNN
	1    1450 3800
	-1   0    0    1   
$EndComp
$Comp
L +15V #PWR4
U 1 1 5AEB69C1
P 1450 3350
F 0 "#PWR4" H 1450 3200 50  0001 C CNN
F 1 "+15V" H 1450 3490 50  0000 C CNN
F 2 "" H 1450 3350 50  0001 C CNN
F 3 "" H 1450 3350 50  0001 C CNN
	1    1450 3350
	1    0    0    -1  
$EndComp
Wire Wire Line
	1450 3350 1450 3450
Wire Wire Line
	1450 3450 1550 3450
Wire Wire Line
	1450 3550 1550 3550
Wire Wire Line
	1550 3650 1450 3650
Wire Wire Line
	1450 3650 1450 3800
Wire Wire Line
	1150 3550 1100 3550
Wire Wire Line
	1100 3550 1100 3600
$Comp
L GND #PWR12
U 1 1 5AEB6FED
P 2450 3950
F 0 "#PWR12" H 2450 3700 50  0001 C CNN
F 1 "GND" H 2450 3800 50  0000 C CNN
F 2 "" H 2450 3950 50  0001 C CNN
F 3 "" H 2450 3950 50  0001 C CNN
	1    2450 3950
	1    0    0    -1  
$EndComp
$Comp
L -15V #PWR11
U 1 1 5AEB6FF3
P 2450 3450
F 0 "#PWR11" H 2450 3550 50  0001 C CNN
F 1 "-15V" H 2450 3600 50  0000 C CNN
F 2 "" H 2450 3450 50  0001 C CNN
F 3 "" H 2450 3450 50  0001 C CNN
	1    2450 3450
	1    0    0    -1  
$EndComp
$Comp
L +15V #PWR15
U 1 1 5AEB6FF9
P 2950 3450
F 0 "#PWR15" H 2950 3300 50  0001 C CNN
F 1 "+15V" H 2950 3590 50  0000 C CNN
F 2 "" H 2950 3450 50  0001 C CNN
F 3 "" H 2950 3450 50  0001 C CNN
	1    2950 3450
	1    0    0    -1  
$EndComp
$Comp
L CP C10
U 1 1 5AEB7242
P 2950 3700
F 0 "C10" H 2975 3800 50  0000 L CNN
F 1 "CP" H 2975 3600 50  0000 L CNN
F 2 "Capacitors_THT:CP_Radial_D5.0mm_P2.00mm" H 2988 3550 50  0001 C CNN
F 3 "" H 2950 3700 50  0001 C CNN
	1    2950 3700
	1    0    0    -1  
$EndComp
$Comp
L CP C7
U 1 1 5AEB72B5
P 2450 3700
F 0 "C7" H 2475 3800 50  0000 L CNN
F 1 "CP" H 2475 3600 50  0000 L CNN
F 2 "Capacitors_THT:CP_Radial_D5.0mm_P2.00mm" H 2488 3550 50  0001 C CNN
F 3 "" H 2450 3700 50  0001 C CNN
	1    2450 3700
	-1   0    0    1   
$EndComp
Wire Wire Line
	2950 3450 2950 3550
Wire Wire Line
	2950 3850 2950 3900
Wire Wire Line
	2950 3900 2450 3900
Wire Wire Line
	2450 3850 2450 3950
Connection ~ 2450 3900
Wire Wire Line
	2450 3550 2450 3450
$Comp
L OP179GS U2
U 1 1 5AEB86F5
P 4350 1600
F 0 "U2" H 4350 1850 50  0000 L CNN
F 1 "OP211" H 4350 1750 50  0000 L CNN
F 2 "Housings_SOIC:SOIC-8_3.9x4.9mm_Pitch1.27mm" H 4350 1600 50  0001 C CNN
F 3 "" H 4500 1750 50  0001 C CNN
	1    4350 1600
	1    0    0    1   
$EndComp
Connection ~ 3350 1500
Wire Wire Line
	3350 1950 3550 1950
Connection ~ 3550 1950
$Comp
L R R9
U 1 1 5AEB9D59
P 3350 1750
F 0 "R9" V 3430 1750 50  0000 C CNN
F 1 "R" V 3350 1750 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 3280 1750 50  0001 C CNN
F 3 "" H 3350 1750 50  0001 C CNN
	1    3350 1750
	-1   0    0    1   
$EndComp
Wire Wire Line
	3350 1600 3350 1500
Wire Wire Line
	3350 1900 3350 1950
$EndSCHEMATC