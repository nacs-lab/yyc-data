EESchema Schematic File Version 2
LIBS:dds_power-rescue
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
LIBS:qd48t018033
LIBS:ptb48510bas
LIBS:dds_power-cache
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
L QD48T018033 U1
U 1 1 556796D6
P 3500 3200
F 0 "U1" H 3500 3150 60  0000 C CNN
F 1 "QD48T018033" H 3500 3450 60  0000 C CNN
F 2 "qd48t018033:qd48t018033" H 3500 3200 60  0001 C CNN
F 3 "" H 3500 3200 60  0000 C CNN
	1    3500 3200
	1    0    0    -1  
$EndComp
$Comp
L R-RESCUE-dds_power R1
U 1 1 55679800
P 2000 3500
F 0 "R1" V 2080 3500 40  0000 C CNN
F 1 "0" V 2007 3501 40  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 1930 3500 30  0001 C CNN
F 3 "" H 2000 3500 30  0000 C CNN
	1    2000 3500
	0    1    1    0   
$EndComp
$Comp
L R-RESCUE-dds_power R2
U 1 1 55679B76
P 4800 3250
F 0 "R2" V 4880 3250 40  0000 C CNN
F 1 "0" V 4807 3251 40  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 4730 3250 30  0001 C CNN
F 3 "" H 4800 3250 30  0000 C CNN
	1    4800 3250
	0    1    1    0   
$EndComp
$Comp
L POT-RESCUE-dds_power RV1
U 1 1 55679BD9
P 4550 2900
F 0 "RV1" H 4550 2800 50  0000 C CNN
F 1 "POT" H 4550 2900 50  0000 C CNN
F 2 "Potentiometers:Potentiometer_Bourns_3296Y_3-8Zoll_Angular_ScrewUp" H 4550 2900 60  0001 C CNN
F 3 "" H 4550 2900 60  0000 C CNN
	1    4550 2900
	0    1    1    0   
$EndComp
$Comp
L CONN_3 K1
U 1 1 5567B186
P 4800 3600
F 0 "K1" V 4750 3600 50  0000 C CNN
F 1 "3.3_1.8" V 4850 3600 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 4800 3600 60  0001 C CNN
F 3 "" H 4800 3600 60  0000 C CNN
	1    4800 3600
	1    0    0    -1  
$EndComp
$Comp
L CONN_3 K2
U 1 1 5567B301
P 4800 4000
F 0 "K2" V 4750 4000 50  0000 C CNN
F 1 "3.3_1.8" V 4850 4000 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 4800 4000 60  0001 C CNN
F 3 "" H 4800 4000 60  0000 C CNN
	1    4800 4000
	1    0    0    -1  
$EndComp
$Comp
L CONN_3 K3
U 1 1 5567B363
P 4800 4400
F 0 "K3" V 4750 4400 50  0000 C CNN
F 1 "3.3_1.8" V 4850 4400 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 4800 4400 60  0001 C CNN
F 3 "" H 4800 4400 60  0000 C CNN
	1    4800 4400
	1    0    0    -1  
$EndComp
$Comp
L CONN_3 K4
U 1 1 5567B3A5
P 4800 4800
F 0 "K4" V 4750 4800 50  0000 C CNN
F 1 "3.3_1.8" V 4850 4800 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 4800 4800 60  0001 C CNN
F 3 "" H 4800 4800 60  0000 C CNN
	1    4800 4800
	1    0    0    -1  
$EndComp
$Comp
L CONN_3 K5
U 1 1 5567B3F0
P 4800 5200
F 0 "K5" V 4750 5200 50  0000 C CNN
F 1 "3.3_1.8" V 4850 5200 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 4800 5200 60  0001 C CNN
F 3 "" H 4800 5200 60  0000 C CNN
	1    4800 5200
	1    0    0    -1  
$EndComp
$Comp
L CONN_3 K6
U 1 1 5567B4C8
P 4800 5600
F 0 "K6" V 4750 5600 50  0000 C CNN
F 1 "3.3_1.8" V 4850 5600 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 4800 5600 60  0001 C CNN
F 3 "" H 4800 5600 60  0000 C CNN
	1    4800 5600
	1    0    0    -1  
$EndComp
$Comp
L C-RESCUE-dds_power C2
U 1 1 55687ADF
P 2050 2800
F 0 "C2" H 2050 2900 40  0000 L CNN
F 1 "C" H 2056 2715 40  0000 L CNN
F 2 "Capacitors_SMD:C_1206" H 2088 2650 30  0001 C CNN
F 3 "" H 2050 2800 60  0000 C CNN
	1    2050 2800
	0    1    1    0   
$EndComp
$Comp
L C-RESCUE-dds_power C4
U 1 1 55687E96
P 4750 6050
F 0 "C4" H 4750 6150 40  0000 L CNN
F 1 "C" H 4756 5965 40  0000 L CNN
F 2 "Capacitors_SMD:C_1206" H 4788 5900 30  0001 C CNN
F 3 "" H 4750 6050 60  0000 C CNN
	1    4750 6050
	1    0    0    -1  
$EndComp
$Comp
L C-RESCUE-dds_power C5
U 1 1 55687EED
P 4750 6550
F 0 "C5" H 4750 6650 40  0000 L CNN
F 1 "C" H 4756 6465 40  0000 L CNN
F 2 "Capacitors_SMD:C_1206" H 4788 6400 30  0001 C CNN
F 3 "" H 4750 6550 60  0000 C CNN
	1    4750 6550
	1    0    0    -1  
$EndComp
$Comp
L CP-RESCUE-dds_power C7
U 1 1 556883ED
P 5000 6550
F 0 "C7" H 5050 6650 40  0000 L CNN
F 1 "CP" H 5050 6450 40  0000 L CNN
F 2 "Capacitors_ThroughHole:C_Radial_D5_L11_P2" H 5100 6400 30  0001 C CNN
F 3 "" H 5000 6550 300 0000 C CNN
	1    5000 6550
	-1   0    0    1   
$EndComp
$Comp
L CP-RESCUE-dds_power C6
U 1 1 556883C0
P 5000 6050
F 0 "C6" H 5050 6150 40  0000 L CNN
F 1 "CP" H 5050 5950 40  0000 L CNN
F 2 "Capacitors_ThroughHole:C_Radial_D5_L11_P2" H 5100 5900 30  0001 C CNN
F 3 "" H 5000 6050 300 0000 C CNN
	1    5000 6050
	1    0    0    -1  
$EndComp
Wire Wire Line
	1600 3050 2800 3050
Wire Wire Line
	1600 3250 1650 3250
Wire Wire Line
	1650 3250 1650 3350
Wire Wire Line
	1650 3350 2800 3350
Wire Wire Line
	1750 3200 2800 3200
Wire Wire Line
	1750 3200 1750 3500
Connection ~ 2250 3350
Wire Wire Line
	4200 3250 4550 3250
Wire Wire Line
	4200 3150 4550 3150
Wire Wire Line
	4700 2900 5050 2900
Wire Wire Line
	5050 2900 5050 3250
Wire Wire Line
	4550 2650 4700 2650
Wire Wire Line
	4700 2650 4700 2900
Wire Wire Line
	7700 3800 7700 3950
Wire Wire Line
	8100 3800 8100 3950
Connection ~ 7700 3800
Wire Wire Line
	8500 3800 8500 3950
Connection ~ 8100 3800
Wire Wire Line
	8900 3800 8900 3950
Connection ~ 8500 3800
Wire Wire Line
	7500 3600 7500 3950
Wire Wire Line
	7900 3600 7900 3950
Connection ~ 7500 3600
Wire Wire Line
	8300 3600 8300 3950
Connection ~ 7900 3600
Wire Wire Line
	8700 3600 8700 3950
Connection ~ 8300 3600
Wire Wire Line
	4200 3050 4400 3050
Wire Wire Line
	4400 3050 4400 6750
Wire Wire Line
	4400 3500 4450 3500
Wire Wire Line
	4400 3900 4450 3900
Connection ~ 4400 3500
Wire Wire Line
	4400 4300 4450 4300
Connection ~ 4400 3900
Wire Wire Line
	4400 4700 4450 4700
Connection ~ 4400 4300
Wire Wire Line
	4400 5100 4450 5100
Connection ~ 4400 4700
Wire Wire Line
	4400 5500 4450 5500
Connection ~ 4400 5100
Wire Wire Line
	4350 3250 4350 6300
Wire Wire Line
	4350 3600 4450 3600
Connection ~ 4350 3250
Wire Wire Line
	4350 4000 4450 4000
Connection ~ 4350 3600
Wire Wire Line
	4350 4400 4450 4400
Connection ~ 4350 4000
Wire Wire Line
	4350 4800 4450 4800
Connection ~ 4350 4400
Wire Wire Line
	4350 5200 4450 5200
Connection ~ 4350 4800
Wire Wire Line
	4350 5600 4450 5600
Connection ~ 4350 5200
Wire Wire Line
	4200 3350 4300 3350
Wire Wire Line
	4300 3350 4300 5850
Wire Wire Line
	4300 3700 4450 3700
Wire Wire Line
	4300 4100 4450 4100
Connection ~ 4300 3700
Wire Wire Line
	4300 4500 4450 4500
Connection ~ 4300 4100
Wire Wire Line
	4300 4900 4450 4900
Connection ~ 4300 4500
Wire Wire Line
	4300 5300 4450 5300
Connection ~ 4300 4900
Wire Wire Line
	4300 5700 4450 5700
Connection ~ 4300 5300
Connection ~ 2250 2800
Connection ~ 1850 3050
Connection ~ 1850 2800
Wire Wire Line
	4300 5850 5000 5850
Connection ~ 4300 5700
Connection ~ 4750 5850
Wire Wire Line
	4750 6250 4750 6350
Wire Wire Line
	5000 6250 5000 6350
Wire Wire Line
	4400 6750 5000 6750
Connection ~ 4400 5500
Connection ~ 4750 6750
Wire Wire Line
	4350 6300 5750 6300
Connection ~ 4750 6300
Connection ~ 5000 6300
Connection ~ 4350 5600
$Comp
L CP-RESCUE-dds_power C1
U 1 1 55688DD5
P 1850 2500
F 0 "C1" H 1900 2600 40  0000 L CNN
F 1 "CP" H 1900 2400 40  0000 L CNN
F 2 "Capacitors_ThroughHole:C_Radial_D5_L11_P2" H 1950 2350 30  0001 C CNN
F 3 "" H 1850 2500 300 0000 C CNN
	1    1850 2500
	-1   0    0    1   
$EndComp
$Comp
L CP-RESCUE-dds_power C3
U 1 1 556890A2
P 2250 2500
F 0 "C3" H 2300 2600 40  0000 L CNN
F 1 "CP" H 2300 2400 40  0000 L CNN
F 2 "Capacitors_ThroughHole:C_Radial_D5_L11_P2" H 2350 2350 30  0001 C CNN
F 3 "" H 2250 2500 300 0000 C CNN
	1    2250 2500
	1    0    0    -1  
$EndComp
Wire Wire Line
	2250 2700 2250 3500
Wire Wire Line
	1850 2700 1850 3050
Wire Wire Line
	2250 2300 1850 2300
$Comp
L CONN_3 P1
U 1 1 55689EDE
P 1250 3150
F 0 "P1" V 1200 3150 50  0000 C CNN
F 1 "CONN_3" V 1300 3150 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 1250 3150 60  0001 C CNN
F 3 "" H 1250 3150 60  0000 C CNN
	1    1250 3150
	-1   0    0    1   
$EndComp
$Comp
L GND-RESCUE-dds_power #PWR01
U 1 1 5568A308
P 1700 3650
F 0 "#PWR01" H 1700 3650 30  0001 C CNN
F 1 "GND" H 1700 3580 30  0001 C CNN
F 2 "" H 1700 3650 60  0000 C CNN
F 3 "" H 1700 3650 60  0000 C CNN
	1    1700 3650
	1    0    0    -1  
$EndComp
Wire Wire Line
	1700 3650 1700 3150
Wire Wire Line
	1700 3150 1600 3150
$Comp
L GND-RESCUE-dds_power #PWR02
U 1 1 5568A57C
P 5750 6400
F 0 "#PWR02" H 5750 6400 30  0001 C CNN
F 1 "GND" H 5750 6330 30  0001 C CNN
F 2 "" H 5750 6400 60  0000 C CNN
F 3 "" H 5750 6400 60  0000 C CNN
	1    5750 6400
	1    0    0    -1  
$EndComp
Wire Wire Line
	5750 6300 5750 6400
$Comp
L PTB48510BAS U2
U 1 1 57ED8DF5
P 3500 2150
F 0 "U2" H 3500 2100 60  0000 C CNN
F 1 "PTB48510BAS" H 3500 2400 60  0000 C CNN
F 2 "ptb48510bas:ptb48510bas" H 3500 2150 60  0001 C CNN
F 3 "" H 3500 2150 60  0000 C CNN
	1    3500 2150
	1    0    0    -1  
$EndComp
NoConn ~ 2800 2100
Wire Wire Line
	2800 2200 2800 2300
Wire Wire Line
	2800 2300 2700 2300
Wire Wire Line
	2700 2300 2700 3350
Connection ~ 2700 3350
Wire Wire Line
	2800 2000 2550 2000
Wire Wire Line
	2550 2000 2550 3050
Connection ~ 2550 3050
$Comp
L POT-RESCUE-dds_power RV2
U 1 1 57EDA2DC
P 4700 2450
F 0 "RV2" H 4700 2350 50  0000 C CNN
F 1 "POT" H 4700 2450 50  0000 C CNN
F 2 "Potentiometers:Potentiometer_Bourns_3296Y_3-8Zoll_Angular_ScrewUp" H 4700 2450 60  0001 C CNN
F 3 "" H 4700 2450 60  0000 C CNN
	1    4700 2450
	1    0    0    -1  
$EndComp
$Comp
L R-RESCUE-dds_power R3
U 1 1 57EDA385
P 4550 2200
F 0 "R3" V 4630 2200 40  0000 C CNN
F 1 "0" V 4557 2201 40  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 4480 2200 30  0001 C CNN
F 3 "" H 4550 2200 30  0000 C CNN
	1    4550 2200
	0    1    1    0   
$EndComp
Wire Wire Line
	4200 2100 5600 2100
Wire Wire Line
	4200 2300 5550 2300
Wire Wire Line
	4200 2200 4300 2200
Wire Wire Line
	4950 2450 4950 2200
Wire Wire Line
	4950 2200 4800 2200
$Comp
L CONN_3 K7
U 1 1 57EDB9DE
P 6050 2550
F 0 "K7" V 6000 2550 50  0000 C CNN
F 1 "+-12" V 6100 2550 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 6050 2550 60  0001 C CNN
F 3 "" H 6050 2550 60  0000 C CNN
	1    6050 2550
	1    0    0    -1  
$EndComp
$Comp
L CONN_3 K8
U 1 1 57EDB9E4
P 6050 2950
F 0 "K8" V 6000 2950 50  0000 C CNN
F 1 "+-12" V 6100 2950 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 6050 2950 60  0001 C CNN
F 3 "" H 6050 2950 60  0000 C CNN
	1    6050 2950
	1    0    0    -1  
$EndComp
$Comp
L CONN_3 K9
U 1 1 57EDB9EA
P 6050 3350
F 0 "K9" V 6000 3350 50  0000 C CNN
F 1 "+-12" V 6100 3350 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 6050 3350 60  0001 C CNN
F 3 "" H 6050 3350 60  0000 C CNN
	1    6050 3350
	1    0    0    -1  
$EndComp
$Comp
L CONN_3 K10
U 1 1 57EDB9F0
P 6050 3750
F 0 "K10" V 6000 3750 50  0000 C CNN
F 1 "+-12" V 6100 3750 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 6050 3750 60  0001 C CNN
F 3 "" H 6050 3750 60  0000 C CNN
	1    6050 3750
	1    0    0    -1  
$EndComp
$Comp
L CONN_3 K11
U 1 1 57EDB9F6
P 6050 4150
F 0 "K11" V 6000 4150 50  0000 C CNN
F 1 "+-12" V 6100 4150 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 6050 4150 60  0001 C CNN
F 3 "" H 6050 4150 60  0000 C CNN
	1    6050 4150
	1    0    0    -1  
$EndComp
$Comp
L CONN_3 K12
U 1 1 57EDB9FC
P 6050 4550
F 0 "K12" V 6000 4550 50  0000 C CNN
F 1 "+-12" V 6100 4550 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 6050 4550 60  0001 C CNN
F 3 "" H 6050 4550 60  0000 C CNN
	1    6050 4550
	1    0    0    -1  
$EndComp
$Comp
L C-RESCUE-dds_power C8
U 1 1 57EDBA02
P 6000 5000
F 0 "C8" H 6000 5100 40  0000 L CNN
F 1 "C" H 6006 4915 40  0000 L CNN
F 2 "Capacitors_SMD:C_1206" H 6038 4850 30  0001 C CNN
F 3 "" H 6000 5000 60  0000 C CNN
	1    6000 5000
	1    0    0    -1  
$EndComp
$Comp
L C-RESCUE-dds_power C9
U 1 1 57EDBA08
P 6000 5500
F 0 "C9" H 6000 5600 40  0000 L CNN
F 1 "C" H 6006 5415 40  0000 L CNN
F 2 "Capacitors_SMD:C_1206" H 6038 5350 30  0001 C CNN
F 3 "" H 6000 5500 60  0000 C CNN
	1    6000 5500
	1    0    0    -1  
$EndComp
$Comp
L CP-RESCUE-dds_power C11
U 1 1 57EDBA0E
P 6250 5500
F 0 "C11" H 6300 5600 40  0000 L CNN
F 1 "CP" H 6300 5400 40  0000 L CNN
F 2 "Capacitors_ThroughHole:C_Radial_D5_L11_P2" H 6350 5350 30  0001 C CNN
F 3 "" H 6250 5500 300 0000 C CNN
	1    6250 5500
	-1   0    0    1   
$EndComp
$Comp
L CP-RESCUE-dds_power C10
U 1 1 57EDBA14
P 6250 5000
F 0 "C10" H 6300 5100 40  0000 L CNN
F 1 "CP" H 6300 4900 40  0000 L CNN
F 2 "Capacitors_ThroughHole:C_Radial_D5_L11_P2" H 6350 4850 30  0001 C CNN
F 3 "" H 6250 5000 300 0000 C CNN
	1    6250 5000
	-1   0    0    1   
$EndComp
Wire Wire Line
	5650 2000 5650 5700
Wire Wire Line
	5650 2450 5700 2450
Wire Wire Line
	5650 2850 5700 2850
Connection ~ 5650 2450
Wire Wire Line
	5650 3250 5700 3250
Connection ~ 5650 2850
Wire Wire Line
	5650 3650 5700 3650
Connection ~ 5650 3250
Wire Wire Line
	5650 4050 5700 4050
Connection ~ 5650 3650
Wire Wire Line
	5650 4450 5700 4450
Connection ~ 5650 4050
Wire Wire Line
	5600 2100 5600 5250
Wire Wire Line
	5600 2550 5700 2550
Wire Wire Line
	5600 2950 5700 2950
Connection ~ 5600 2550
Wire Wire Line
	5600 3350 5700 3350
Connection ~ 5600 2950
Wire Wire Line
	5600 3750 5700 3750
Connection ~ 5600 3350
Wire Wire Line
	5600 4150 5700 4150
Connection ~ 5600 3750
Wire Wire Line
	5600 4550 5700 4550
Connection ~ 5600 4150
Wire Wire Line
	5550 2300 5550 4800
Wire Wire Line
	5550 2650 5700 2650
Wire Wire Line
	5550 3050 5700 3050
Connection ~ 5550 2650
Wire Wire Line
	5550 3450 5700 3450
Connection ~ 5550 3050
Wire Wire Line
	5550 3850 5700 3850
Connection ~ 5550 3450
Wire Wire Line
	5550 4250 5700 4250
Connection ~ 5550 3850
Wire Wire Line
	5550 4650 5700 4650
Connection ~ 5550 4250
Wire Wire Line
	5550 4800 6250 4800
Connection ~ 5550 4650
Connection ~ 6000 4800
Wire Wire Line
	6000 5200 6000 5300
Wire Wire Line
	6250 5200 6250 5300
Wire Wire Line
	5650 5700 6250 5700
Connection ~ 5650 4450
Connection ~ 6000 5700
Wire Wire Line
	5600 5250 7000 5250
Connection ~ 6000 5250
Connection ~ 6250 5250
Connection ~ 5600 4550
$Comp
L GND-RESCUE-dds_power #PWR03
U 1 1 57EDBA53
P 7000 5350
F 0 "#PWR03" H 7000 5350 30  0001 C CNN
F 1 "GND" H 7000 5280 30  0001 C CNN
F 2 "" H 7000 5350 60  0000 C CNN
F 3 "" H 7000 5350 60  0000 C CNN
	1    7000 5350
	1    0    0    -1  
$EndComp
Wire Wire Line
	7000 5250 7000 5350
Wire Wire Line
	4200 2000 5650 2000
Connection ~ 4700 2300
Wire Wire Line
	4450 2450 4450 2300
Connection ~ 4450 2300
$Comp
L CONN_3 K14
U 1 1 57EDDBFB
P 7000 3700
F 0 "K14" V 6950 3700 50  0000 C CNN
F 1 "+-5" V 7050 3700 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 7000 3700 60  0001 C CNN
F 3 "" H 7000 3700 60  0000 C CNN
	1    7000 3700
	-1   0    0    1   
$EndComp
$Comp
L CONN_3 K16
U 1 1 57EDE18B
P 7600 4300
F 0 "K16" V 7550 4300 50  0000 C CNN
F 1 "+-5" V 7650 4300 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 7600 4300 60  0001 C CNN
F 3 "" H 7600 4300 60  0000 C CNN
	1    7600 4300
	0    1    1    0   
$EndComp
$Comp
L CONN_3 K18
U 1 1 57EDE253
P 8000 4300
F 0 "K18" V 7950 4300 50  0000 C CNN
F 1 "+-5" V 8050 4300 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 8000 4300 60  0001 C CNN
F 3 "" H 8000 4300 60  0000 C CNN
	1    8000 4300
	0    1    1    0   
$EndComp
$Comp
L CONN_3 K20
U 1 1 57EDE2E5
P 8400 4300
F 0 "K20" V 8350 4300 50  0000 C CNN
F 1 "+-5" V 8450 4300 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 8400 4300 60  0001 C CNN
F 3 "" H 8400 4300 60  0000 C CNN
	1    8400 4300
	0    1    1    0   
$EndComp
$Comp
L CONN_3 K22
U 1 1 57EDE37A
P 8800 4300
F 0 "K22" V 8750 4300 50  0000 C CNN
F 1 "+-5" V 8850 4300 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 8800 4300 60  0001 C CNN
F 3 "" H 8800 4300 60  0000 C CNN
	1    8800 4300
	0    1    1    0   
$EndComp
Wire Wire Line
	7600 3950 7600 3700
Connection ~ 7600 3700
Wire Wire Line
	8000 3700 8000 3950
Connection ~ 8000 3700
Wire Wire Line
	8400 3700 8400 3950
Connection ~ 8400 3700
Wire Wire Line
	8800 3700 8800 3950
Wire Wire Line
	10300 3600 7350 3600
Wire Wire Line
	7350 3700 9800 3700
Wire Wire Line
	9300 3800 7350 3800
$Comp
L C-RESCUE-dds_power C15
U 1 1 57EE0AF9
P 9500 3300
F 0 "C15" H 9500 3400 40  0000 L CNN
F 1 "C" H 9506 3215 40  0000 L CNN
F 2 "Capacitors_SMD:C_1206" H 9538 3150 30  0001 C CNN
F 3 "" H 9500 3300 60  0000 C CNN
	1    9500 3300
	0    1    1    0   
$EndComp
$Comp
L C-RESCUE-dds_power C19
U 1 1 57EE0BF0
P 10100 3300
F 0 "C19" H 10100 3400 40  0000 L CNN
F 1 "C" H 10106 3215 40  0000 L CNN
F 2 "Capacitors_SMD:C_1206" H 10138 3150 30  0001 C CNN
F 3 "" H 10100 3300 60  0000 C CNN
	1    10100 3300
	0    1    1    0   
$EndComp
$Comp
L CP-RESCUE-dds_power C14
U 1 1 57EE0CB8
P 9500 3000
F 0 "C14" H 9550 3100 40  0000 L CNN
F 1 "CP" H 9550 2900 40  0000 L CNN
F 2 "Capacitors_ThroughHole:C_Radial_D5_L11_P2" H 9600 2850 30  0001 C CNN
F 3 "" H 9500 3000 300 0000 C CNN
	1    9500 3000
	0    -1   -1   0   
$EndComp
$Comp
L CP-RESCUE-dds_power C18
U 1 1 57EE1403
P 10100 3000
F 0 "C18" H 10150 3100 40  0000 L CNN
F 1 "CP" H 10150 2900 40  0000 L CNN
F 2 "Capacitors_ThroughHole:C_Radial_D5_L11_P2" H 10200 2850 30  0001 C CNN
F 3 "" H 10100 3000 300 0000 C CNN
	1    10100 3000
	0    -1   -1   0   
$EndComp
Wire Wire Line
	9800 3000 9800 3950
Wire Wire Line
	9700 3000 9900 3000
Connection ~ 8800 3700
Connection ~ 9800 3000
Wire Wire Line
	9700 3300 9900 3300
Connection ~ 9800 3300
Wire Wire Line
	9300 3000 9300 3800
Connection ~ 8900 3800
Connection ~ 9300 3300
Wire Wire Line
	10300 3000 10300 3600
Connection ~ 8700 3600
Connection ~ 10300 3300
$Comp
L GND-RESCUE-dds_power #PWR04
U 1 1 57EE1DEE
P 9800 3950
F 0 "#PWR04" H 9800 3950 30  0001 C CNN
F 1 "GND" H 9800 3880 30  0001 C CNN
F 2 "" H 9800 3950 60  0000 C CNN
F 3 "" H 9800 3950 60  0000 C CNN
	1    9800 3950
	1    0    0    -1  
$EndComp
Connection ~ 9800 3700
Wire Wire Line
	7700 2200 7700 2350
Wire Wire Line
	8100 2200 8100 2350
Connection ~ 7700 2200
Wire Wire Line
	8500 2200 8500 2350
Connection ~ 8100 2200
Wire Wire Line
	8900 2200 8900 2350
Connection ~ 8500 2200
Wire Wire Line
	7500 2000 7500 2350
Wire Wire Line
	7900 2000 7900 2350
Connection ~ 7500 2000
Wire Wire Line
	8300 2000 8300 2350
Connection ~ 7900 2000
Wire Wire Line
	8700 2000 8700 2350
Connection ~ 8300 2000
$Comp
L CONN_3 K13
U 1 1 57EE2578
P 7000 2100
F 0 "K13" V 6950 2100 50  0000 C CNN
F 1 "+-5" V 7050 2100 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 7000 2100 60  0001 C CNN
F 3 "" H 7000 2100 60  0000 C CNN
	1    7000 2100
	-1   0    0    1   
$EndComp
$Comp
L CONN_3 K15
U 1 1 57EE257E
P 7600 2700
F 0 "K15" V 7550 2700 50  0000 C CNN
F 1 "+-5" V 7650 2700 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 7600 2700 60  0001 C CNN
F 3 "" H 7600 2700 60  0000 C CNN
	1    7600 2700
	0    1    1    0   
$EndComp
$Comp
L CONN_3 K17
U 1 1 57EE2584
P 8000 2700
F 0 "K17" V 7950 2700 50  0000 C CNN
F 1 "+-5" V 8050 2700 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 8000 2700 60  0001 C CNN
F 3 "" H 8000 2700 60  0000 C CNN
	1    8000 2700
	0    1    1    0   
$EndComp
$Comp
L CONN_3 K19
U 1 1 57EE258A
P 8400 2700
F 0 "K19" V 8350 2700 50  0000 C CNN
F 1 "+-5" V 8450 2700 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 8400 2700 60  0001 C CNN
F 3 "" H 8400 2700 60  0000 C CNN
	1    8400 2700
	0    1    1    0   
$EndComp
$Comp
L CONN_3 K21
U 1 1 57EE2590
P 8800 2700
F 0 "K21" V 8750 2700 50  0000 C CNN
F 1 "+-5" V 8850 2700 40  0000 C CNN
F 2 "cmr_link:3-pins_150mil" H 8800 2700 60  0001 C CNN
F 3 "" H 8800 2700 60  0000 C CNN
	1    8800 2700
	0    1    1    0   
$EndComp
Wire Wire Line
	7600 2350 7600 2100
Connection ~ 7600 2100
Wire Wire Line
	8000 2100 8000 2350
Connection ~ 8000 2100
Wire Wire Line
	8400 2100 8400 2350
Connection ~ 8400 2100
Wire Wire Line
	8800 2100 8800 2350
Wire Wire Line
	10300 2000 7350 2000
Wire Wire Line
	7350 2100 9800 2100
Wire Wire Line
	9300 2200 7350 2200
$Comp
L C-RESCUE-dds_power C13
U 1 1 57EE25A0
P 9500 1700
F 0 "C13" H 9500 1800 40  0000 L CNN
F 1 "C" H 9506 1615 40  0000 L CNN
F 2 "Capacitors_SMD:C_1206" H 9538 1550 30  0001 C CNN
F 3 "" H 9500 1700 60  0000 C CNN
	1    9500 1700
	0    1    1    0   
$EndComp
$Comp
L C-RESCUE-dds_power C17
U 1 1 57EE25A6
P 10100 1700
F 0 "C17" H 10100 1800 40  0000 L CNN
F 1 "C" H 10106 1615 40  0000 L CNN
F 2 "Capacitors_SMD:C_1206" H 10138 1550 30  0001 C CNN
F 3 "" H 10100 1700 60  0000 C CNN
	1    10100 1700
	0    1    1    0   
$EndComp
$Comp
L CP-RESCUE-dds_power C12
U 1 1 57EE25AC
P 9500 1400
F 0 "C12" H 9550 1500 40  0000 L CNN
F 1 "CP" H 9550 1300 40  0000 L CNN
F 2 "Capacitors_ThroughHole:C_Radial_D5_L11_P2" H 9600 1250 30  0001 C CNN
F 3 "" H 9500 1400 300 0000 C CNN
	1    9500 1400
	0    -1   -1   0   
$EndComp
$Comp
L CP-RESCUE-dds_power C16
U 1 1 57EE25B2
P 10100 1400
F 0 "C16" H 10150 1500 40  0000 L CNN
F 1 "CP" H 10150 1300 40  0000 L CNN
F 2 "Capacitors_ThroughHole:C_Radial_D5_L11_P2" H 10200 1250 30  0001 C CNN
F 3 "" H 10100 1400 300 0000 C CNN
	1    10100 1400
	0    -1   -1   0   
$EndComp
Wire Wire Line
	9800 1400 9800 2350
Wire Wire Line
	9700 1400 9900 1400
Connection ~ 8800 2100
Connection ~ 9800 1400
Wire Wire Line
	9700 1700 9900 1700
Connection ~ 9800 1700
Wire Wire Line
	9300 1400 9300 2200
Connection ~ 8900 2200
Connection ~ 9300 1700
Wire Wire Line
	10300 1400 10300 2000
Connection ~ 8700 2000
Connection ~ 10300 1700
$Comp
L GND-RESCUE-dds_power #PWR05
U 1 1 57EE25C4
P 9800 2350
F 0 "#PWR05" H 9800 2350 30  0001 C CNN
F 1 "GND" H 9800 2280 30  0001 C CNN
F 2 "" H 9800 2350 60  0000 C CNN
F 3 "" H 9800 2350 60  0000 C CNN
	1    9800 2350
	1    0    0    -1  
$EndComp
Connection ~ 9800 2100
$EndSCHEMATC
