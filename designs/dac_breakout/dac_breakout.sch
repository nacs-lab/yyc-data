EESchema Schematic File Version 2
LIBS:power
LIBS:device
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
LIBS:rj45x8
LIBS:dac_breakout-cache
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
L CONN_02X30 P?
U 1 1 58A77384
P 2500 3200
F 0 "P?" H 2500 4750 50  0000 C CNN
F 1 "CONN_02X30" V 2500 3200 50  0000 C CNN
F 2 "" H 2500 2700 50  0000 C CNN
F 3 "" H 2500 2700 50  0000 C CNN
	1    2500 3200
	1    0    0    -1  
$EndComp
$Comp
L RJ45x8 U?
U 1 1 58AB30FD
P 5550 3550
F 0 "U?" H 5550 3500 60  0000 C CNN
F 1 "RJ45x8" H 5550 3600 60  0000 C CNN
F 2 "" H 5550 3550 60  0001 C CNN
F 3 "" H 5550 3550 60  0001 C CNN
	1    5550 3550
	0    -1   -1   0   
$EndComp
$Comp
L GND #PWR?
U 1 1 58AB32AD
P 5250 6350
F 0 "#PWR?" H 5250 6100 50  0001 C CNN
F 1 "GND" H 5250 6200 50  0000 C CNN
F 2 "" H 5250 6350 50  0000 C CNN
F 3 "" H 5250 6350 50  0000 C CNN
	1    5250 6350
	1    0    0    -1  
$EndComp
Wire Wire Line
	5250 6250 5850 6250
Connection ~ 5350 6250
Connection ~ 5450 6250
Connection ~ 5550 6250
Connection ~ 5650 6250
Connection ~ 5750 6250
Wire Wire Line
	5250 6250 5250 6350
$Comp
L CONN_01X03 P?
U 1 1 58AB3423
P 2000 5700
F 0 "P?" H 2000 5900 50  0000 C CNN
F 1 "CONN_01X03" V 2100 5700 50  0000 C CNN
F 2 "" H 2000 5700 50  0000 C CNN
F 3 "" H 2000 5700 50  0000 C CNN
	1    2000 5700
	1    0    0    -1  
$EndComp
$Comp
L CONN_01X03 P?
U 1 1 58AB3483
P 2000 6200
F 0 "P?" H 2000 6400 50  0000 C CNN
F 1 "CONN_01X03" V 2100 6200 50  0000 C CNN
F 2 "" H 2000 6200 50  0000 C CNN
F 3 "" H 2000 6200 50  0000 C CNN
	1    2000 6200
	1    0    0    -1  
$EndComp
$EndSCHEMATC
