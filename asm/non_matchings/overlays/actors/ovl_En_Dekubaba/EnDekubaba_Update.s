glabel EnDekubaba_Update
/* 02C20 809E83F0 27BDFFC8 */  addiu   $sp, $sp, 0xFFC8           ## $sp = FFFFFFC8
/* 02C24 809E83F4 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 02C28 809E83F8 AFB10020 */  sw      $s1, 0x0020($sp)           
/* 02C2C 809E83FC AFB0001C */  sw      $s0, 0x001C($sp)           
/* 02C30 809E8400 90820248 */  lbu     $v0, 0x0248($a0)           ## 00000248
/* 02C34 809E8404 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 02C38 809E8408 00A08825 */  or      $s1, $a1, $zero            ## $s1 = 00000000
/* 02C3C 809E840C 304E0002 */  andi    $t6, $v0, 0x0002           ## $t6 = 00000000
/* 02C40 809E8410 11C00003 */  beq     $t6, $zero, .L809E8420     
/* 02C44 809E8414 304FFFFD */  andi    $t7, $v0, 0xFFFD           ## $t7 = 00000000
/* 02C48 809E8418 0C27981E */  jal     func_809E6078              
/* 02C4C 809E841C A08F0248 */  sb      $t7, 0x0248($a0)           ## 00000248
.L809E8420:
/* 02C50 809E8420 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 02C54 809E8424 0C27A050 */  jal     func_809E8140              
/* 02C58 809E8428 02202825 */  or      $a1, $s1, $zero            ## $a1 = 00000000
/* 02C5C 809E842C 8E1901C0 */  lw      $t9, 0x01C0($s0)           ## 000001C0
/* 02C60 809E8430 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 02C64 809E8434 02202825 */  or      $a1, $s1, $zero            ## $a1 = 00000000
/* 02C68 809E8438 0320F809 */  jalr    $ra, $t9                   
/* 02C6C 809E843C 00000000 */  nop
/* 02C70 809E8440 8E0201C0 */  lw      $v0, 0x01C0($s0)           ## 000001C0
/* 02C74 809E8444 3C18809E */  lui     $t8, %hi(func_809E7BB0)    ## $t8 = 809E0000
/* 02C78 809E8448 27187BB0 */  addiu   $t8, $t8, %lo(func_809E7BB0) ## $t8 = 809E7BB0
/* 02C7C 809E844C 17020015 */  bne     $t8, $v0, .L809E84A4       
/* 02C80 809E8450 3C03809F */  lui     $v1, %hi(func_809E80D8)    ## $v1 = 809F0000
/* 02C84 809E8454 0C00B638 */  jal     Actor_MoveForward
              
/* 02C88 809E8458 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 02C8C 809E845C 3C014120 */  lui     $at, 0x4120                ## $at = 41200000
/* 02C90 809E8460 44810000 */  mtc1    $at, $f0                   ## $f0 = 10.00
/* 02C94 809E8464 3C014170 */  lui     $at, 0x4170                ## $at = 41700000
/* 02C98 809E8468 44813000 */  mtc1    $at, $f6                   ## $f6 = 15.00
/* 02C9C 809E846C C6040230 */  lwc1    $f4, 0x0230($s0)           ## 00000230
/* 02CA0 809E8470 24080005 */  addiu   $t0, $zero, 0x0005         ## $t0 = 00000005
/* 02CA4 809E8474 44060000 */  mfc1    $a2, $f0                   
/* 02CA8 809E8478 46062202 */  mul.s   $f8, $f4, $f6              
/* 02CAC 809E847C AFA80014 */  sw      $t0, 0x0014($sp)           
/* 02CB0 809E8480 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 02CB4 809E8484 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 02CB8 809E8488 E7A00010 */  swc1    $f0, 0x0010($sp)           
/* 02CBC 809E848C 44074000 */  mfc1    $a3, $f8                   
/* 02CC0 809E8490 0C00B92D */  jal     func_8002E4B4              
/* 02CC4 809E8494 00000000 */  nop
/* 02CC8 809E8498 3C03809F */  lui     $v1, %hi(func_809E80D8)    ## $v1 = 809F0000
/* 02CCC 809E849C 10000013 */  beq     $zero, $zero, .L809E84EC   
/* 02CD0 809E84A0 246380D8 */  addiu   $v1, $v1, %lo(func_809E80D8) ## $v1 = 809E80D8
.L809E84A4:
/* 02CD4 809E84A4 246380D8 */  addiu   $v1, $v1, %lo(func_809E80D8) ## $v1 = 809E01B0
/* 02CD8 809E84A8 10620010 */  beq     $v1, $v0, .L809E84EC       
/* 02CDC 809E84AC 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 02CE0 809E84B0 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 02CE4 809E84B4 24090004 */  addiu   $t1, $zero, 0x0004         ## $t1 = 00000004
/* 02CE8 809E84B8 AFA90014 */  sw      $t1, 0x0014($sp)           
/* 02CEC 809E84BC 44060000 */  mfc1    $a2, $f0                   
/* 02CF0 809E84C0 44070000 */  mfc1    $a3, $f0                   
/* 02CF4 809E84C4 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 02CF8 809E84C8 AFA3002C */  sw      $v1, 0x002C($sp)           
/* 02CFC 809E84CC 0C00B92D */  jal     func_8002E4B4              
/* 02D00 809E84D0 E7A00010 */  swc1    $f0, 0x0010($sp)           
/* 02D04 809E84D4 8E0A0234 */  lw      $t2, 0x0234($s0)           ## 00000234
/* 02D08 809E84D8 8FA3002C */  lw      $v1, 0x002C($sp)           
/* 02D0C 809E84DC 55400004 */  bnel    $t2, $zero, .L809E84F0     
/* 02D10 809E84E0 8E0D01C0 */  lw      $t5, 0x01C0($s0)           ## 000001C0
/* 02D14 809E84E4 8E0B0078 */  lw      $t3, 0x0078($s0)           ## 00000078
/* 02D18 809E84E8 AE0B0234 */  sw      $t3, 0x0234($s0)           ## 00000234
.L809E84EC:
/* 02D1C 809E84EC 8E0D01C0 */  lw      $t5, 0x01C0($s0)           ## 000001C0
.L809E84F0:
/* 02D20 809E84F0 3C0C809E */  lui     $t4, %hi(func_809E7104)    ## $t4 = 809E0000
/* 02D24 809E84F4 258C7104 */  addiu   $t4, $t4, %lo(func_809E7104) ## $t4 = 809E7104
/* 02D28 809E84F8 158D000C */  bne     $t4, $t5, .L809E852C       
/* 02D2C 809E84FC 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 02D30 809E8500 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 02D34 809E8504 34211E60 */  ori     $at, $at, 0x1E60           ## $at = 00011E60
/* 02D38 809E8508 02212821 */  addu    $a1, $s1, $at              
/* 02D3C 809E850C 26060238 */  addiu   $a2, $s0, 0x0238           ## $a2 = 00000238
/* 02D40 809E8510 0C0175E7 */  jal     Actor_CollisionCheck_SetAT
              ## CollisionCheck_setAT
/* 02D44 809E8514 AFA3002C */  sw      $v1, 0x002C($sp)           
/* 02D48 809E8518 8E0E0004 */  lw      $t6, 0x0004($s0)           ## 00000004
/* 02D4C 809E851C 3C010100 */  lui     $at, 0x0100                ## $at = 01000000
/* 02D50 809E8520 8FA3002C */  lw      $v1, 0x002C($sp)           
/* 02D54 809E8524 01C17825 */  or      $t7, $t6, $at              ## $t7 = 01000000
/* 02D58 809E8528 AE0F0004 */  sw      $t7, 0x0004($s0)           ## 00000004
.L809E852C:
/* 02D5C 809E852C 92190249 */  lbu     $t9, 0x0249($s0)           ## 00000249
/* 02D60 809E8530 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 02D64 809E8534 34211E60 */  ori     $at, $at, 0x1E60           ## $at = 00011E60
/* 02D68 809E8538 33380001 */  andi    $t8, $t9, 0x0001           ## $t8 = 00000000
/* 02D6C 809E853C 13000006 */  beq     $t8, $zero, .L809E8558     
/* 02D70 809E8540 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 02D74 809E8544 02212821 */  addu    $a1, $s1, $at              
/* 02D78 809E8548 26060238 */  addiu   $a2, $s0, 0x0238           ## $a2 = 00000238
/* 02D7C 809E854C 0C01767D */  jal     Actor_CollisionCheck_SetAC
              ## CollisionCheck_setAC
/* 02D80 809E8550 AFA3002C */  sw      $v1, 0x002C($sp)           
/* 02D84 809E8554 8FA3002C */  lw      $v1, 0x002C($sp)           
.L809E8558:
/* 02D88 809E8558 8E0801C0 */  lw      $t0, 0x01C0($s0)           ## 000001C0
/* 02D8C 809E855C 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 02D90 809E8560 34211E60 */  ori     $at, $at, 0x1E60           ## $at = 00011E60
/* 02D94 809E8564 10680004 */  beq     $v1, $t0, .L809E8578       
/* 02D98 809E8568 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 02D9C 809E856C 02212821 */  addu    $a1, $s1, $at              
/* 02DA0 809E8570 0C017713 */  jal     Actor_CollisionCheck_SetOT
              ## CollisionCheck_setOT
/* 02DA4 809E8574 26060238 */  addiu   $a2, $s0, 0x0238           ## $a2 = 00000238
.L809E8578:
/* 02DA8 809E8578 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 02DAC 809E857C 8FB0001C */  lw      $s0, 0x001C($sp)           
/* 02DB0 809E8580 8FB10020 */  lw      $s1, 0x0020($sp)           
/* 02DB4 809E8584 03E00008 */  jr      $ra                        
/* 02DB8 809E8588 27BD0038 */  addiu   $sp, $sp, 0x0038           ## $sp = 00000000


