glabel func_80A33444
/* 00874 80A33444 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 00878 80A33448 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 0087C 80A3344C AFB00020 */  sw      $s0, 0x0020($sp)           
/* 00880 80A33450 948202F4 */  lhu     $v0, 0x02F4($a0)           ## 000002F4
/* 00884 80A33454 44802000 */  mtc1    $zero, $f4                 ## $f4 = 0.00
/* 00888 80A33458 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 0088C 80A3345C 304E0010 */  andi    $t6, $v0, 0x0010           ## $t6 = 00000000
/* 00890 80A33460 00A03025 */  or      $a2, $a1, $zero            ## $a2 = 00000000
/* 00894 80A33464 11C00004 */  beq     $t6, $zero, .L80A33478     
/* 00898 80A33468 E4840068 */  swc1    $f4, 0x0068($a0)           ## 00000068
/* 0089C 80A3346C 304FFFEF */  andi    $t7, $v0, 0xFFEF           ## $t7 = 00000000
/* 008A0 80A33470 10000018 */  beq     $zero, $zero, .L80A334D4   
/* 008A4 80A33474 A48F02F4 */  sh      $t7, 0x02F4($a0)           ## 000002F4
.L80A33478:
/* 008A8 80A33478 00C02025 */  or      $a0, $a2, $zero            ## $a0 = 00000000
/* 008AC 80A3347C 0C28CBB3 */  jal     func_80A32ECC              
/* 008B0 80A33480 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 008B4 80A33484 10400009 */  beq     $v0, $zero, .L80A334AC     
/* 008B8 80A33488 24180064 */  addiu   $t8, $zero, 0x0064         ## $t8 = 00000064
/* 008BC 80A3348C 92080306 */  lbu     $t0, 0x0306($s0)           ## 00000306
/* 008C0 80A33490 8619008A */  lh      $t9, 0x008A($s0)           ## 0000008A
/* 008C4 80A33494 A2180305 */  sb      $t8, 0x0305($s0)           ## 00000305
/* 008C8 80A33498 0102082A */  slt     $at, $t0, $v0              
/* 008CC 80A3349C 1020000D */  beq     $at, $zero, .L80A334D4     
/* 008D0 80A334A0 A61902F8 */  sh      $t9, 0x02F8($s0)           ## 000002F8
/* 008D4 80A334A4 1000000B */  beq     $zero, $zero, .L80A334D4   
/* 008D8 80A334A8 A2020306 */  sb      $v0, 0x0306($s0)           ## 00000306
.L80A334AC:
/* 008DC 80A334AC 860902F8 */  lh      $t1, 0x02F8($s0)           ## 000002F8
/* 008E0 80A334B0 860A0032 */  lh      $t2, 0x0032($s0)           ## 00000032
/* 008E4 80A334B4 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 008E8 80A334B8 24050001 */  addiu   $a1, $zero, 0x0001         ## $a1 = 00000001
/* 008EC 80A334BC 552A0006 */  bnel    $t1, $t2, .L80A334D8       
/* 008F0 80A334C0 92020306 */  lbu     $v0, 0x0306($s0)           ## 00000306
/* 008F4 80A334C4 0C28CAF4 */  jal     func_80A32BD0              
/* 008F8 80A334C8 A2000306 */  sb      $zero, 0x0306($s0)         ## 00000306
/* 008FC 80A334CC 1000001B */  beq     $zero, $zero, .L80A3353C   
/* 00900 80A334D0 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80A334D4:
/* 00904 80A334D4 92020306 */  lbu     $v0, 0x0306($s0)           ## 00000306
.L80A334D8:
/* 00908 80A334D8 24010001 */  addiu   $at, $zero, 0x0001         ## $at = 00000001
/* 0090C 80A334DC 26040032 */  addiu   $a0, $s0, 0x0032           ## $a0 = 00000032
/* 00910 80A334E0 10410006 */  beq     $v0, $at, .L80A334FC       
/* 00914 80A334E4 24060002 */  addiu   $a2, $zero, 0x0002         ## $a2 = 00000002
/* 00918 80A334E8 24010002 */  addiu   $at, $zero, 0x0002         ## $at = 00000002
/* 0091C 80A334EC 1041000A */  beq     $v0, $at, .L80A33518       
/* 00920 80A334F0 26040032 */  addiu   $a0, $s0, 0x0032           ## $a0 = 00000032
/* 00924 80A334F4 1000000F */  beq     $zero, $zero, .L80A33534   
/* 00928 80A334F8 86030032 */  lh      $v1, 0x0032($s0)           ## 00000032
.L80A334FC:
/* 0092C 80A334FC 860502F8 */  lh      $a1, 0x02F8($s0)           ## 000002F8
/* 00930 80A33500 240B0100 */  addiu   $t3, $zero, 0x0100         ## $t3 = 00000100
/* 00934 80A33504 AFAB0010 */  sw      $t3, 0x0010($sp)           
/* 00938 80A33508 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS
              
/* 0093C 80A3350C 24070200 */  addiu   $a3, $zero, 0x0200         ## $a3 = 00000200
/* 00940 80A33510 10000008 */  beq     $zero, $zero, .L80A33534   
/* 00944 80A33514 86030032 */  lh      $v1, 0x0032($s0)           ## 00000032
.L80A33518:
/* 00948 80A33518 860502F8 */  lh      $a1, 0x02F8($s0)           ## 000002F8
/* 0094C 80A3351C 240C0180 */  addiu   $t4, $zero, 0x0180         ## $t4 = 00000180
/* 00950 80A33520 AFAC0010 */  sw      $t4, 0x0010($sp)           
/* 00954 80A33524 24060002 */  addiu   $a2, $zero, 0x0002         ## $a2 = 00000002
/* 00958 80A33528 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS
              
/* 0095C 80A3352C 24070600 */  addiu   $a3, $zero, 0x0600         ## $a3 = 00000600
/* 00960 80A33530 86030032 */  lh      $v1, 0x0032($s0)           ## 00000032
.L80A33534:
/* 00964 80A33534 A60300B6 */  sh      $v1, 0x00B6($s0)           ## 000000B6
/* 00968 80A33538 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80A3353C:
/* 0096C 80A3353C 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 00970 80A33540 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 00974 80A33544 03E00008 */  jr      $ra                        
/* 00978 80A33548 00000000 */  nop


