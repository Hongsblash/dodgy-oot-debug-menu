glabel func_809E18B0
/* 00000 809E18B0 27BDFFC8 */  addiu   $sp, $sp, 0xFFC8           ## $sp = FFFFFFC8
/* 00004 809E18B4 AFBF0024 */  sw      $ra, 0x0024($sp)
/* 00008 809E18B8 AFA40038 */  sw      $a0, 0x0038($sp)
/* 0000C 809E18BC AFA60040 */  sw      $a2, 0x0040($sp)
/* 00010 809E18C0 8CC20000 */  lw      $v0, 0x0000($a2)           ## 00000000
/* 00014 809E18C4 44803000 */  mtc1    $zero, $f6                 ## $f6 = 0.00
/* 00018 809E18C8 00057900 */  sll     $t7, $a1,  4
/* 0001C 809E18CC 04400003 */  bltz    $v0, .L809E18DC
/* 00020 809E18D0 3C18809E */  lui     $t8, %hi(D_809E2918)       ## $t8 = 809E0000
/* 00024 809E18D4 14A20005 */  bne     $a1, $v0, .L809E18EC
/* 00028 809E18D8 0005C900 */  sll     $t9, $a1,  4
.L809E18DC:
/* 0002C 809E18DC 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 00030 809E18E0 27182918 */  addiu   $t8, $t8, %lo(D_809E2918)  ## $t8 = 809E2918
/* 00034 809E18E4 10000005 */  beq     $zero, $zero, .L809E18FC
/* 00038 809E18E8 01F81821 */  addu    $v1, $t7, $t8
.L809E18EC:
/* 0003C 809E18EC 3C08809E */  lui     $t0, %hi(D_809E2918)       ## $t0 = 809E0000
/* 00040 809E18F0 25082918 */  addiu   $t0, $t0, %lo(D_809E2918)  ## $t0 = 809E2918
/* 00044 809E18F4 03281821 */  addu    $v1, $t9, $t0
/* 00048 809E18F8 C460000C */  lwc1    $f0, 0x000C($v1)           ## 0000000C
.L809E18FC:
/* 0004C 809E18FC C4640004 */  lwc1    $f4, 0x0004($v1)           ## 00000004
/* 00050 809E1900 4604303E */  c.le.s  $f6, $f4
/* 00054 809E1904 00000000 */  nop
/* 00058 809E1908 45020017 */  bc1fl   .L809E1968
/* 0005C 809E190C 8C640000 */  lw      $a0, 0x0000($v1)           ## 00000000
/* 00060 809E1910 8C640000 */  lw      $a0, 0x0000($v1)           ## 00000000
/* 00064 809E1914 AFA3002C */  sw      $v1, 0x002C($sp)
/* 00068 809E1918 AFA5003C */  sw      $a1, 0x003C($sp)
/* 0006C 809E191C 0C028800 */  jal     SkelAnime_GetFrameCount

/* 00070 809E1920 E7A00034 */  swc1    $f0, 0x0034($sp)
/* 00074 809E1924 44824000 */  mtc1    $v0, $f8                   ## $f8 = 0.00
/* 00078 809E1928 8FA3002C */  lw      $v1, 0x002C($sp)
/* 0007C 809E192C C7A00034 */  lwc1    $f0, 0x0034($sp)
/* 00080 809E1930 468042A0 */  cvt.s.w $f10, $f8
/* 00084 809E1934 8FA40038 */  lw      $a0, 0x0038($sp)
/* 00088 809E1938 90690008 */  lbu     $t1, 0x0008($v1)           ## 00000008
/* 0008C 809E193C 24070000 */  addiu   $a3, $zero, 0x0000         ## $a3 = 00000000
/* 00090 809E1940 8C650000 */  lw      $a1, 0x0000($v1)           ## 00000000
/* 00094 809E1944 8C660004 */  lw      $a2, 0x0004($v1)           ## 00000004
/* 00098 809E1948 E7AA0010 */  swc1    $f10, 0x0010($sp)
/* 0009C 809E194C E7A00018 */  swc1    $f0, 0x0018($sp)
/* 000A0 809E1950 2484014C */  addiu   $a0, $a0, 0x014C           ## $a0 = 0000014C
/* 000A4 809E1954 0C029468 */  jal     SkelAnime_ChangeAnim

/* 000A8 809E1958 AFA90014 */  sw      $t1, 0x0014($sp)
/* 000AC 809E195C 10000016 */  beq     $zero, $zero, .L809E19B8
/* 000B0 809E1960 8FAB003C */  lw      $t3, 0x003C($sp)
/* 000B4 809E1964 8C640000 */  lw      $a0, 0x0000($v1)           ## 00000000
.L809E1968:
/* 000B8 809E1968 AFA3002C */  sw      $v1, 0x002C($sp)
/* 000BC 809E196C AFA5003C */  sw      $a1, 0x003C($sp)
/* 000C0 809E1970 0C028800 */  jal     SkelAnime_GetFrameCount

/* 000C4 809E1974 E7A00034 */  swc1    $f0, 0x0034($sp)
/* 000C8 809E1978 44828000 */  mtc1    $v0, $f16                  ## $f16 = 0.00
/* 000CC 809E197C 8FA3002C */  lw      $v1, 0x002C($sp)
/* 000D0 809E1980 C7A00034 */  lwc1    $f0, 0x0034($sp)
/* 000D4 809E1984 46808420 */  cvt.s.w $f16, $f16
/* 000D8 809E1988 8FA40038 */  lw      $a0, 0x0038($sp)
/* 000DC 809E198C 44809000 */  mtc1    $zero, $f18                ## $f18 = 0.00
/* 000E0 809E1990 906A0008 */  lbu     $t2, 0x0008($v1)           ## 00000008
/* 000E4 809E1994 8C650000 */  lw      $a1, 0x0000($v1)           ## 00000000
/* 000E8 809E1998 8C660004 */  lw      $a2, 0x0004($v1)           ## 00000004
/* 000EC 809E199C 44078000 */  mfc1    $a3, $f16
/* 000F0 809E19A0 E7A00018 */  swc1    $f0, 0x0018($sp)
/* 000F4 809E19A4 2484014C */  addiu   $a0, $a0, 0x014C           ## $a0 = 0000014C
/* 000F8 809E19A8 AFAA0014 */  sw      $t2, 0x0014($sp)
/* 000FC 809E19AC 0C029468 */  jal     SkelAnime_ChangeAnim

/* 00100 809E19B0 E7B20010 */  swc1    $f18, 0x0010($sp)
/* 00104 809E19B4 8FAB003C */  lw      $t3, 0x003C($sp)
.L809E19B8:
/* 00108 809E19B8 8FAC0040 */  lw      $t4, 0x0040($sp)
/* 0010C 809E19BC AD8B0000 */  sw      $t3, 0x0000($t4)           ## 00000000
/* 00110 809E19C0 8FBF0024 */  lw      $ra, 0x0024($sp)
/* 00114 809E19C4 27BD0038 */  addiu   $sp, $sp, 0x0038           ## $sp = 00000000
/* 00118 809E19C8 03E00008 */  jr      $ra
/* 0011C 809E19CC 00000000 */  nop


