glabel func_80A3F300
/* 00D90 80A3F300 27BDFFC8 */  addiu   $sp, $sp, 0xFFC8           ## $sp = FFFFFFC8
/* 00D94 80A3F304 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 00D98 80A3F308 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 00D9C 80A3F30C AFA5003C */  sw      $a1, 0x003C($sp)           
/* 00DA0 80A3F310 8486001C */  lh      $a2, 0x001C($a0)           ## 0000001C
/* 00DA4 80A3F314 2401000F */  addiu   $at, $zero, 0x000F         ## $at = 0000000F
/* 00DA8 80A3F318 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 00DAC 80A3F31C 30C6000F */  andi    $a2, $a2, 0x000F           ## $a2 = 00000000
/* 00DB0 80A3F320 14C10003 */  bne     $a2, $at, .L80A3F330       
/* 00DB4 80A3F324 8FAE003C */  lw      $t6, 0x003C($sp)           
/* 00DB8 80A3F328 1000005F */  beq     $zero, $zero, .L80A3F4A8   
/* 00DBC 80A3F32C 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
.L80A3F330:
/* 00DC0 80A3F330 3C0F0001 */  lui     $t7, 0x0001                ## $t7 = 00010000
/* 00DC4 80A3F334 01EE7821 */  addu    $t7, $t7, $t6              
/* 00DC8 80A3F338 8DEF1E08 */  lw      $t7, 0x1E08($t7)           ## 00011E08
/* 00DCC 80A3F33C 0006C0C0 */  sll     $t8, $a2,  3               
/* 00DD0 80A3F340 3C0A8016 */  lui     $t2, %hi(gSegments)
/* 00DD4 80A3F344 01F82021 */  addu    $a0, $t7, $t8              
/* 00DD8 80A3F348 8C830004 */  lw      $v1, 0x0004($a0)           ## 00000004
/* 00DDC 80A3F34C 860C0218 */  lh      $t4, 0x0218($s0)           ## 00000218
/* 00DE0 80A3F350 3C0100FF */  lui     $at, 0x00FF                ## $at = 00FF0000
/* 00DE4 80A3F354 0003C900 */  sll     $t9, $v1,  4               
/* 00DE8 80A3F358 00194702 */  srl     $t0, $t9, 28               
/* 00DEC 80A3F35C 00084880 */  sll     $t1, $t0,  2               
/* 00DF0 80A3F360 01495021 */  addu    $t2, $t2, $t1              
/* 00DF4 80A3F364 8D4A6FA8 */  lw      $t2, %lo(gSegments)($t2)
/* 00DF8 80A3F368 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = 00FFFFFF
/* 00DFC 80A3F36C 00615824 */  and     $t3, $v1, $at              
/* 00E00 80A3F370 000C6880 */  sll     $t5, $t4,  2               
/* 00E04 80A3F374 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 00E08 80A3F378 01AC6823 */  subu    $t5, $t5, $t4              
/* 00E0C 80A3F37C 014B1021 */  addu    $v0, $t2, $t3              
/* 00E10 80A3F380 00411021 */  addu    $v0, $v0, $at              
/* 00E14 80A3F384 000D6840 */  sll     $t5, $t5,  1               
/* 00E18 80A3F388 004D1021 */  addu    $v0, $v0, $t5              
/* 00E1C 80A3F38C 844E0000 */  lh      $t6, 0x0000($v0)           ## 00000000
/* 00E20 80A3F390 844F0004 */  lh      $t7, 0x0004($v0)           ## 00000004
/* 00E24 80A3F394 C6080024 */  lwc1    $f8, 0x0024($s0)           ## 00000024
/* 00E28 80A3F398 448E2000 */  mtc1    $t6, $f4                   ## $f4 = 0.00
/* 00E2C 80A3F39C 448F5000 */  mtc1    $t7, $f10                  ## $f10 = 0.00
/* 00E30 80A3F3A0 C612002C */  lwc1    $f18, 0x002C($s0)          ## 0000002C
/* 00E34 80A3F3A4 468021A0 */  cvt.s.w $f6, $f4                   
/* 00E38 80A3F3A8 AFA40034 */  sw      $a0, 0x0034($sp)           
/* 00E3C 80A3F3AC 46805420 */  cvt.s.w $f16, $f10                 
/* 00E40 80A3F3B0 46083301 */  sub.s   $f12, $f6, $f8             
/* 00E44 80A3F3B4 46128381 */  sub.s   $f14, $f16, $f18           
/* 00E48 80A3F3B8 E7AC002C */  swc1    $f12, 0x002C($sp)          
/* 00E4C 80A3F3BC 0C03F494 */  jal     Math_FAtan2F              
/* 00E50 80A3F3C0 E7AE0028 */  swc1    $f14, 0x0028($sp)          
/* 00E54 80A3F3C4 3C0180A4 */  lui     $at, %hi(D_80A420D8)       ## $at = 80A40000
/* 00E58 80A3F3C8 C42420D8 */  lwc1    $f4, %lo(D_80A420D8)($at)  
/* 00E5C 80A3F3CC 24190001 */  addiu   $t9, $zero, 0x0001         ## $t9 = 00000001
/* 00E60 80A3F3D0 AFB90010 */  sw      $t9, 0x0010($sp)           
/* 00E64 80A3F3D4 46040182 */  mul.s   $f6, $f0, $f4              
/* 00E68 80A3F3D8 26040032 */  addiu   $a0, $s0, 0x0032           ## $a0 = 00000032
/* 00E6C 80A3F3DC 2406000A */  addiu   $a2, $zero, 0x000A         ## $a2 = 0000000A
/* 00E70 80A3F3E0 240703E8 */  addiu   $a3, $zero, 0x03E8         ## $a3 = 000003E8
/* 00E74 80A3F3E4 4600320D */  trunc.w.s $f8, $f6                   
/* 00E78 80A3F3E8 44054000 */  mfc1    $a1, $f8                   
/* 00E7C 80A3F3EC 00000000 */  nop
/* 00E80 80A3F3F0 00052C00 */  sll     $a1, $a1, 16               
/* 00E84 80A3F3F4 0C01E1A7 */  jal     Math_SmoothStepToS
              
/* 00E88 80A3F3F8 00052C03 */  sra     $a1, $a1, 16               
/* 00E8C 80A3F3FC C7AC002C */  lwc1    $f12, 0x002C($sp)          
/* 00E90 80A3F400 C7AE0028 */  lwc1    $f14, 0x0028($sp)          
/* 00E94 80A3F404 3C014416 */  lui     $at, 0x4416                ## $at = 44160000
/* 00E98 80A3F408 460C6282 */  mul.s   $f10, $f12, $f12           
/* 00E9C 80A3F40C 44812000 */  mtc1    $at, $f4                   ## $f4 = 600.00
/* 00EA0 80A3F410 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 00EA4 80A3F414 460E7402 */  mul.s   $f16, $f14, $f14           
/* 00EA8 80A3F418 46105480 */  add.s   $f18, $f10, $f16           
/* 00EAC 80A3F41C 4604903C */  c.lt.s  $f18, $f4                  
/* 00EB0 80A3F420 00000000 */  nop
/* 00EB4 80A3F424 45000020 */  bc1f    .L80A3F4A8                 
/* 00EB8 80A3F428 00000000 */  nop
/* 00EBC 80A3F42C 86080218 */  lh      $t0, 0x0218($s0)           ## 00000218
/* 00EC0 80A3F430 25090001 */  addiu   $t1, $t0, 0x0001           ## $t1 = 00000001
/* 00EC4 80A3F434 A6090218 */  sh      $t1, 0x0218($s0)           ## 00000218
/* 00EC8 80A3F438 8FAB0034 */  lw      $t3, 0x0034($sp)           
/* 00ECC 80A3F43C 860A0218 */  lh      $t2, 0x0218($s0)           ## 00000218
/* 00ED0 80A3F440 916C0000 */  lbu     $t4, 0x0000($t3)           ## 00000000
/* 00ED4 80A3F444 014C082A */  slt     $at, $t2, $t4              
/* 00ED8 80A3F448 54200003 */  bnel    $at, $zero, .L80A3F458     
/* 00EDC 80A3F44C 8602001C */  lh      $v0, 0x001C($s0)           ## 0000001C
/* 00EE0 80A3F450 A6000218 */  sh      $zero, 0x0218($s0)         ## 00000218
/* 00EE4 80A3F454 8602001C */  lh      $v0, 0x001C($s0)           ## 0000001C
.L80A3F458:
/* 00EE8 80A3F458 8FA4003C */  lw      $a0, 0x003C($sp)           
/* 00EEC 80A3F45C 304D00F0 */  andi    $t5, $v0, 0x00F0           ## $t5 = 00000000
/* 00EF0 80A3F460 11A00003 */  beq     $t5, $zero, .L80A3F470     
/* 00EF4 80A3F464 00000000 */  nop
/* 00EF8 80A3F468 1000000F */  beq     $zero, $zero, .L80A3F4A8   
/* 00EFC 80A3F46C 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
.L80A3F470:
/* 00F00 80A3F470 0C00B2D0 */  jal     Flags_GetSwitch
              
/* 00F04 80A3F474 00022A03 */  sra     $a1, $v0,  8               
/* 00F08 80A3F478 50400004 */  beql    $v0, $zero, .L80A3F48C     
/* 00F0C 80A3F47C 860E0218 */  lh      $t6, 0x0218($s0)           ## 00000218
/* 00F10 80A3F480 10000009 */  beq     $zero, $zero, .L80A3F4A8   
/* 00F14 80A3F484 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
/* 00F18 80A3F488 860E0218 */  lh      $t6, 0x0218($s0)           ## 00000218
.L80A3F48C:
/* 00F1C 80A3F48C 860F00B8 */  lh      $t7, 0x00B8($s0)           ## 000000B8
/* 00F20 80A3F490 01CF082A */  slt     $at, $t6, $t7              
/* 00F24 80A3F494 14200002 */  bne     $at, $zero, .L80A3F4A0     
/* 00F28 80A3F498 00000000 */  nop
/* 00F2C 80A3F49C A6000218 */  sh      $zero, 0x0218($s0)         ## 00000218
.L80A3F4A0:
/* 00F30 80A3F4A0 10000001 */  beq     $zero, $zero, .L80A3F4A8   
/* 00F34 80A3F4A4 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
.L80A3F4A8:
/* 00F38 80A3F4A8 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 00F3C 80A3F4AC 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 00F40 80A3F4B0 27BD0038 */  addiu   $sp, $sp, 0x0038           ## $sp = 00000000
/* 00F44 80A3F4B4 03E00008 */  jr      $ra                        
/* 00F48 80A3F4B8 00000000 */  nop
