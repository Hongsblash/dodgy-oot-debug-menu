glabel func_80996D14
/* 00A74 80996D14 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 00A78 80996D18 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 00A7C 80996D1C AFB00018 */  sw      $s0, 0x0018($sp)           
/* 00A80 80996D20 AFA50024 */  sw      $a1, 0x0024($sp)           
/* 00A84 80996D24 908E016C */  lbu     $t6, 0x016C($a0)           ## 0000016C
/* 00A88 80996D28 24010003 */  addiu   $at, $zero, 0x0003         ## $at = 00000003
/* 00A8C 80996D2C 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 00A90 80996D30 51C1001E */  beql    $t6, $at, .L80996DAC       
/* 00A94 80996D34 860F0166 */  lh      $t7, 0x0166($s0)           ## 00000166
/* 00A98 80996D38 44802000 */  mtc1    $zero, $f4                 ## $f4 = 0.00
/* 00A9C 80996D3C C4860060 */  lwc1    $f6, 0x0060($a0)           ## 00000060
/* 00AA0 80996D40 46062032 */  c.eq.s  $f4, $f6                   
/* 00AA4 80996D44 00000000 */  nop
/* 00AA8 80996D48 45020007 */  bc1fl   .L80996D68                 
/* 00AAC 80996D4C 26040060 */  addiu   $a0, $s0, 0x0060           ## $a0 = 00000060
/* 00AB0 80996D50 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 00AB4 80996D54 24052814 */  addiu   $a1, $zero, 0x2814         ## $a1 = 00002814
/* 00AB8 80996D58 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00ABC 80996D5C 0C265B18 */  jal     func_80996C60              
/* 00AC0 80996D60 8FA50024 */  lw      $a1, 0x0024($sp)           
/* 00AC4 80996D64 26040060 */  addiu   $a0, $s0, 0x0060           ## $a0 = 00000060
.L80996D68:
/* 00AC8 80996D68 3C054170 */  lui     $a1, 0x4170                ## $a1 = 41700000
/* 00ACC 80996D6C 0C01DE80 */  jal     Math_ApproxF
              
/* 00AD0 80996D70 3C064040 */  lui     $a2, 0x4040                ## $a2 = 40400000
/* 00AD4 80996D74 3C014348 */  lui     $at, 0x4348                ## $at = 43480000
/* 00AD8 80996D78 44815000 */  mtc1    $at, $f10                  ## $f10 = 200.00
/* 00ADC 80996D7C C608000C */  lwc1    $f8, 0x000C($s0)           ## 0000000C
/* 00AE0 80996D80 26040028 */  addiu   $a0, $s0, 0x0028           ## $a0 = 00000028
/* 00AE4 80996D84 8E060060 */  lw      $a2, 0x0060($s0)           ## 00000060
/* 00AE8 80996D88 460A4400 */  add.s   $f16, $f8, $f10            
/* 00AEC 80996D8C 44058000 */  mfc1    $a1, $f16                  
/* 00AF0 80996D90 0C01DE80 */  jal     Math_ApproxF
              
/* 00AF4 80996D94 00000000 */  nop
/* 00AF8 80996D98 50400016 */  beql    $v0, $zero, .L80996DF4     
/* 00AFC 80996D9C 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 00B00 80996DA0 10000014 */  beq     $zero, $zero, .L80996DF4   
/* 00B04 80996DA4 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
/* 00B08 80996DA8 860F0166 */  lh      $t7, 0x0166($s0)           ## 00000166
.L80996DAC:
/* 00B0C 80996DAC 24010064 */  addiu   $at, $zero, 0x0064         ## $at = 00000064
/* 00B10 80996DB0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00B14 80996DB4 55E10007 */  bnel    $t7, $at, .L80996DD4       
/* 00B18 80996DB8 26040166 */  addiu   $a0, $s0, 0x0166           ## $a0 = 00000166
/* 00B1C 80996DBC 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 00B20 80996DC0 24052864 */  addiu   $a1, $zero, 0x2864         ## $a1 = 00002864
/* 00B24 80996DC4 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00B28 80996DC8 0C265B18 */  jal     func_80996C60              
/* 00B2C 80996DCC 8FA50024 */  lw      $a1, 0x0024($sp)           
/* 00B30 80996DD0 26040166 */  addiu   $a0, $s0, 0x0166           ## $a0 = 00000166
.L80996DD4:
/* 00B34 80996DD4 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 00B38 80996DD8 0C01DE5F */  jal     Math_ApproxS
              
/* 00B3C 80996DDC 2406000A */  addiu   $a2, $zero, 0x000A         ## $a2 = 0000000A
/* 00B40 80996DE0 50400004 */  beql    $v0, $zero, .L80996DF4     
/* 00B44 80996DE4 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 00B48 80996DE8 10000002 */  beq     $zero, $zero, .L80996DF4   
/* 00B4C 80996DEC 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
/* 00B50 80996DF0 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
.L80996DF4:
/* 00B54 80996DF4 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 00B58 80996DF8 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 00B5C 80996DFC 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 00B60 80996E00 03E00008 */  jr      $ra                        
/* 00B64 80996E04 00000000 */  nop


