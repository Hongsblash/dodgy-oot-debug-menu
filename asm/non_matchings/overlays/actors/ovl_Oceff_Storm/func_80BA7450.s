glabel func_80BA7450
/* 00370 80BA7450 27BDFF80 */  addiu   $sp, $sp, 0xFF80           ## $sp = FFFFFF80
/* 00374 80BA7454 AFBF0034 */  sw      $ra, 0x0034($sp)           
/* 00378 80BA7458 AFA40080 */  sw      $a0, 0x0080($sp)           
/* 0037C 80BA745C AFA50084 */  sw      $a1, 0x0084($sp)           
/* 00380 80BA7460 8CAF009C */  lw      $t7, 0x009C($a1)           ## 0000009C
/* 00384 80BA7464 3C0680BB */  lui     $a2, %hi(D_80BA8C20)       ## $a2 = 80BB0000
/* 00388 80BA7468 24C68C20 */  addiu   $a2, $a2, %lo(D_80BA8C20)  ## $a2 = 80BA8C20
/* 0038C 80BA746C 31F80FFF */  andi    $t8, $t7, 0x0FFF           ## $t8 = 00000000
/* 00390 80BA7470 AFB8007C */  sw      $t8, 0x007C($sp)           
/* 00394 80BA7474 8CA50000 */  lw      $a1, 0x0000($a1)           ## 00000000
/* 00398 80BA7478 27A40064 */  addiu   $a0, $sp, 0x0064           ## $a0 = FFFFFFE4
/* 0039C 80BA747C 240701C1 */  addiu   $a3, $zero, 0x01C1         ## $a3 = 000001C1
/* 003A0 80BA7480 0C031AB1 */  jal     func_800C6AC4              
/* 003A4 80BA7484 AFA50074 */  sw      $a1, 0x0074($sp)           
/* 003A8 80BA7488 8FA30074 */  lw      $v1, 0x0074($sp)           
/* 003AC 80BA748C 3C0BE700 */  lui     $t3, 0xE700                ## $t3 = E7000000
/* 003B0 80BA7490 8C6202D0 */  lw      $v0, 0x02D0($v1)           ## 000002D0
/* 003B4 80BA7494 24590008 */  addiu   $t9, $v0, 0x0008           ## $t9 = 00000008
/* 003B8 80BA7498 AC7902D0 */  sw      $t9, 0x02D0($v1)           ## 000002D0
/* 003BC 80BA749C AC400004 */  sw      $zero, 0x0004($v0)         ## 00000004
/* 003C0 80BA74A0 AC4B0000 */  sw      $t3, 0x0000($v0)           ## 00000000
/* 003C4 80BA74A4 8C6402D0 */  lw      $a0, 0x02D0($v1)           ## 000002D0
/* 003C8 80BA74A8 0C024FCD */  jal     func_80093F34              
/* 003CC 80BA74AC AFA30074 */  sw      $v1, 0x0074($sp)           
/* 003D0 80BA74B0 8FA30074 */  lw      $v1, 0x0074($sp)           
/* 003D4 80BA74B4 244C0008 */  addiu   $t4, $v0, 0x0008           ## $t4 = 00000008
/* 003D8 80BA74B8 3C0DE300 */  lui     $t5, 0xE300                ## $t5 = E3000000
/* 003DC 80BA74BC AC6202D0 */  sw      $v0, 0x02D0($v1)           ## 000002D0
/* 003E0 80BA74C0 AC6C02D0 */  sw      $t4, 0x02D0($v1)           ## 000002D0
/* 003E4 80BA74C4 35AD1A01 */  ori     $t5, $t5, 0x1A01           ## $t5 = E3001A01
/* 003E8 80BA74C8 240F0020 */  addiu   $t7, $zero, 0x0020         ## $t7 = 00000020
/* 003EC 80BA74CC AC4F0004 */  sw      $t7, 0x0004($v0)           ## 00000004
/* 003F0 80BA74D0 AC4D0000 */  sw      $t5, 0x0000($v0)           ## 00000000
/* 003F4 80BA74D4 8C6202D0 */  lw      $v0, 0x02D0($v1)           ## 000002D0
/* 003F8 80BA74D8 3C0EE300 */  lui     $t6, 0xE300                ## $t6 = E3000000
/* 003FC 80BA74DC 35CE1801 */  ori     $t6, $t6, 0x1801           ## $t6 = E3001801
/* 00400 80BA74E0 24580008 */  addiu   $t8, $v0, 0x0008           ## $t8 = 00000008
/* 00404 80BA74E4 AC7802D0 */  sw      $t8, 0x02D0($v1)           ## 000002D0
/* 00408 80BA74E8 24190080 */  addiu   $t9, $zero, 0x0080         ## $t9 = 00000080
/* 0040C 80BA74EC AC590004 */  sw      $t9, 0x0004($v0)           ## 00000004
/* 00410 80BA74F0 AC4E0000 */  sw      $t6, 0x0000($v0)           ## 00000000
/* 00414 80BA74F4 8C6202D0 */  lw      $v0, 0x02D0($v1)           ## 000002D0
/* 00418 80BA74F8 3C0CFA00 */  lui     $t4, 0xFA00                ## $t4 = FA000000
/* 0041C 80BA74FC 358C8080 */  ori     $t4, $t4, 0x8080           ## $t4 = FA008080
/* 00420 80BA7500 244B0008 */  addiu   $t3, $v0, 0x0008           ## $t3 = 00000008
/* 00424 80BA7504 AC6B02D0 */  sw      $t3, 0x02D0($v1)           ## 000002D0
/* 00428 80BA7508 AC4C0000 */  sw      $t4, 0x0000($v0)           ## 00000000
/* 0042C 80BA750C 8FAD0080 */  lw      $t5, 0x0080($sp)           
/* 00430 80BA7510 3C01C8C8 */  lui     $at, 0xC8C8                ## $at = C8C80000
/* 00434 80BA7514 34219600 */  ori     $at, $at, 0x9600           ## $at = C8C89600
/* 00438 80BA7518 91B8014E */  lbu     $t8, 0x014E($t5)           ## E3001B4F
/* 0043C 80BA751C 3C0B80BB */  lui     $t3, %hi(D_80BA8890)       ## $t3 = 80BB0000
/* 00440 80BA7520 256B8890 */  addiu   $t3, $t3, %lo(D_80BA8890)  ## $t3 = 80BA8890
/* 00444 80BA7524 03017025 */  or      $t6, $t8, $at              ## $t6 = C8C89608
/* 00448 80BA7528 AC4E0004 */  sw      $t6, 0x0004($v0)           ## 00000004
/* 0044C 80BA752C 8C6202D0 */  lw      $v0, 0x02D0($v1)           ## 000002D0
/* 00450 80BA7530 3C09DE00 */  lui     $t1, 0xDE00                ## $t1 = DE000000
/* 00454 80BA7534 240E0001 */  addiu   $t6, $zero, 0x0001         ## $t6 = 00000001
/* 00458 80BA7538 24590008 */  addiu   $t9, $v0, 0x0008           ## $t9 = 00000008
/* 0045C 80BA753C AC7902D0 */  sw      $t9, 0x02D0($v1)           ## 000002D0
/* 00460 80BA7540 AC490000 */  sw      $t1, 0x0000($v0)           ## 00000000
/* 00464 80BA7544 AC4B0004 */  sw      $t3, 0x0004($v0)           ## 00000004
/* 00468 80BA7548 8C6202D0 */  lw      $v0, 0x02D0($v1)           ## 000002D0
/* 0046C 80BA754C 8FAA007C */  lw      $t2, 0x007C($sp)           
/* 00470 80BA7550 240B0040 */  addiu   $t3, $zero, 0x0040         ## $t3 = 00000040
/* 00474 80BA7554 244C0008 */  addiu   $t4, $v0, 0x0008           ## $t4 = 00000008
/* 00478 80BA7558 AC6C02D0 */  sw      $t4, 0x02D0($v1)           ## 000002D0
/* 0047C 80BA755C AC490000 */  sw      $t1, 0x0000($v0)           ## 00000000
/* 00480 80BA7560 8FAD0084 */  lw      $t5, 0x0084($sp)           
/* 00484 80BA7564 24190040 */  addiu   $t9, $zero, 0x0040         ## $t9 = 00000040
/* 00488 80BA7568 24180040 */  addiu   $t8, $zero, 0x0040         ## $t8 = 00000040
/* 0048C 80BA756C 8DA40000 */  lw      $a0, 0x0000($t5)           ## E3001A01
/* 00490 80BA7570 240F0040 */  addiu   $t7, $zero, 0x0040         ## $t7 = 00000040
/* 00494 80BA7574 000A3880 */  sll     $a3, $t2,  2               
/* 00498 80BA7578 AFA7001C */  sw      $a3, 0x001C($sp)           
/* 0049C 80BA757C AFA70020 */  sw      $a3, 0x0020($sp)           
/* 004A0 80BA7580 AFAF0010 */  sw      $t7, 0x0010($sp)           
/* 004A4 80BA7584 AFAB0028 */  sw      $t3, 0x0028($sp)           
/* 004A8 80BA7588 AFB90024 */  sw      $t9, 0x0024($sp)           
/* 004AC 80BA758C AFAE0018 */  sw      $t6, 0x0018($sp)           
/* 004B0 80BA7590 AFB80014 */  sw      $t8, 0x0014($sp)           
/* 004B4 80BA7594 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 004B8 80BA7598 AFA30074 */  sw      $v1, 0x0074($sp)           
/* 004BC 80BA759C AFA2004C */  sw      $v0, 0x004C($sp)           
/* 004C0 80BA75A0 0C0253D0 */  jal     Draw_TwoTexScroll              
/* 004C4 80BA75A4 000A30C0 */  sll     $a2, $t2,  3               
/* 004C8 80BA75A8 8FA8004C */  lw      $t0, 0x004C($sp)           
/* 004CC 80BA75AC 8FA30074 */  lw      $v1, 0x0074($sp)           
/* 004D0 80BA75B0 3C0DE450 */  lui     $t5, 0xE450                ## $t5 = E4500000
/* 004D4 80BA75B4 AD020004 */  sw      $v0, 0x0004($t0)           ## 00000004
/* 004D8 80BA75B8 8C6202D0 */  lw      $v0, 0x02D0($v1)           ## 000002D0
/* 004DC 80BA75BC 35AD03C0 */  ori     $t5, $t5, 0x03C0           ## $t5 = E45003C0
/* 004E0 80BA75C0 3C18E100 */  lui     $t8, 0xE100                ## $t8 = E1000000
/* 004E4 80BA75C4 244C0008 */  addiu   $t4, $v0, 0x0008           ## $t4 = 00000008
/* 004E8 80BA75C8 AC6C02D0 */  sw      $t4, 0x02D0($v1)           ## 000002D0
/* 004EC 80BA75CC AC400004 */  sw      $zero, 0x0004($v0)         ## 00000004
/* 004F0 80BA75D0 AC4D0000 */  sw      $t5, 0x0000($v0)           ## 00000000
/* 004F4 80BA75D4 8C6202D0 */  lw      $v0, 0x02D0($v1)           ## 000002D0
/* 004F8 80BA75D8 3C0B008C */  lui     $t3, 0x008C                ## $t3 = 008C0000
/* 004FC 80BA75DC 356BFF74 */  ori     $t3, $t3, 0xFF74           ## $t3 = 008CFF74
/* 00500 80BA75E0 244F0008 */  addiu   $t7, $v0, 0x0008           ## $t7 = 00000008
/* 00504 80BA75E4 AC6F02D0 */  sw      $t7, 0x02D0($v1)           ## 000002D0
/* 00508 80BA75E8 AC400004 */  sw      $zero, 0x0004($v0)         ## 00000004
/* 0050C 80BA75EC AC580000 */  sw      $t8, 0x0000($v0)           ## 00000000
/* 00510 80BA75F0 8C6202D0 */  lw      $v0, 0x02D0($v1)           ## 000002D0
/* 00514 80BA75F4 3C19F100 */  lui     $t9, 0xF100                ## $t9 = F1000000
/* 00518 80BA75F8 3C0680BB */  lui     $a2, %hi(D_80BA8C34)       ## $a2 = 80BB0000
/* 0051C 80BA75FC 244E0008 */  addiu   $t6, $v0, 0x0008           ## $t6 = 00000008
/* 00520 80BA7600 AC6E02D0 */  sw      $t6, 0x02D0($v1)           ## 000002D0
/* 00524 80BA7604 AC4B0004 */  sw      $t3, 0x0004($v0)           ## 00000004
/* 00528 80BA7608 AC590000 */  sw      $t9, 0x0000($v0)           ## 00000000
/* 0052C 80BA760C 8FAC0084 */  lw      $t4, 0x0084($sp)           
/* 00530 80BA7610 24C68C34 */  addiu   $a2, $a2, %lo(D_80BA8C34)  ## $a2 = 80BA8C34
/* 00534 80BA7614 27A40064 */  addiu   $a0, $sp, 0x0064           ## $a0 = FFFFFFE4
/* 00538 80BA7618 240701DD */  addiu   $a3, $zero, 0x01DD         ## $a3 = 000001DD
/* 0053C 80BA761C 0C031AD5 */  jal     func_800C6B54              
/* 00540 80BA7620 8D850000 */  lw      $a1, 0x0000($t4)           ## 00000008
/* 00544 80BA7624 8FBF0034 */  lw      $ra, 0x0034($sp)           
/* 00548 80BA7628 27BD0080 */  addiu   $sp, $sp, 0x0080           ## $sp = 00000000
/* 0054C 80BA762C 03E00008 */  jr      $ra                        
/* 00550 80BA7630 00000000 */  nop


