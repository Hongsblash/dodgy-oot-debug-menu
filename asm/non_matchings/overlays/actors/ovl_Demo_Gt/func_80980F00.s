glabel func_80980F00
/* 03890 80980F00 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 03894 80980F04 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 03898 80980F08 3C014120 */  lui     $at, 0x4120                ## $at = 41200000
/* 0389C 80980F0C 44810000 */  mtc1    $at, $f0                   ## $f0 = 10.00
/* 038A0 80980F10 C4840050 */  lwc1    $f4, 0x0050($a0)           ## 00000050
/* 038A4 80980F14 C4880054 */  lwc1    $f8, 0x0054($a0)           ## 00000054
/* 038A8 80980F18 C4900058 */  lwc1    $f16, 0x0058($a0)          ## 00000058
/* 038AC 80980F1C 46002182 */  mul.s   $f6, $f4, $f0              
/* 038B0 80980F20 24060003 */  addiu   $a2, $zero, 0x0003         ## $a2 = 00000003
/* 038B4 80980F24 24070004 */  addiu   $a3, $zero, 0x0004         ## $a3 = 00000004
/* 038B8 80980F28 46004282 */  mul.s   $f10, $f8, $f0             
/* 038BC 80980F2C 00000000 */  nop
/* 038C0 80980F30 46008482 */  mul.s   $f18, $f16, $f0            
/* 038C4 80980F34 E4860050 */  swc1    $f6, 0x0050($a0)           ## 00000050
/* 038C8 80980F38 E48A0054 */  swc1    $f10, 0x0054($a0)          ## 00000054
/* 038CC 80980F3C E4920058 */  swc1    $f18, 0x0058($a0)          ## 00000058
/* 038D0 80980F40 0C25FB91 */  jal     func_8097EE44              
/* 038D4 80980F44 AFA00010 */  sw      $zero, 0x0010($sp)         
/* 038D8 80980F48 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 038DC 80980F4C 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 038E0 80980F50 03E00008 */  jr      $ra                        
/* 038E4 80980F54 00000000 */  nop


