.rdata
glabel D_80B116A0
    .asciz "\n\n"
    .balign 4

glabel D_80B116A4
    .asciz "[32m☆☆☆☆☆ 親父登場！！むほほほほほほほーん ☆☆☆☆☆ \n[m"
    .balign 4

.text
glabel EnSyatekiMan_Init
/* 00000 80B10870 27BDFFD0 */  addiu   $sp, $sp, 0xFFD0           ## $sp = FFFFFFD0
/* 00004 80B10874 AFB00028 */  sw      $s0, 0x0028($sp)
/* 00008 80B10878 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 0000C 80B1087C AFBF002C */  sw      $ra, 0x002C($sp)
/* 00010 80B10880 3C0480B1 */  lui     $a0, %hi(D_80B116A0)       ## $a0 = 80B10000
/* 00014 80B10884 AFA50034 */  sw      $a1, 0x0034($sp)
/* 00018 80B10888 0C00084C */  jal     osSyncPrintf

/* 0001C 80B1088C 248416A0 */  addiu   $a0, $a0, %lo(D_80B116A0)  ## $a0 = 80B116A0
/* 00020 80B10890 3C0480B1 */  lui     $a0, %hi(D_80B116A4)       ## $a0 = 80B10000
/* 00024 80B10894 0C00084C */  jal     osSyncPrintf

/* 00028 80B10898 248416A4 */  addiu   $a0, $a0, %lo(D_80B116A4)  ## $a0 = 80B116A4
/* 0002C 80B1089C 240E0001 */  addiu   $t6, $zero, 0x0001         ## $t6 = 00000001
/* 00030 80B108A0 3C053C23 */  lui     $a1, 0x3C23                ## $a1 = 3C230000
/* 00034 80B108A4 A20E001F */  sb      $t6, 0x001F($s0)           ## 0000001F
/* 00038 80B108A8 34A5D70A */  ori     $a1, $a1, 0xD70A           ## $a1 = 3C23D70A
/* 0003C 80B108AC 0C00B58B */  jal     Actor_SetScale

/* 00040 80B108B0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00044 80B108B4 3C060601 */  lui     $a2, 0x0601                ## $a2 = 06010000
/* 00048 80B108B8 3C070600 */  lui     $a3, 0x0600                ## $a3 = 06000000
/* 0004C 80B108BC 260F0190 */  addiu   $t7, $s0, 0x0190           ## $t7 = 00000190
/* 00050 80B108C0 261801C6 */  addiu   $t8, $s0, 0x01C6           ## $t8 = 000001C6
/* 00054 80B108C4 24190009 */  addiu   $t9, $zero, 0x0009         ## $t9 = 00000009
/* 00058 80B108C8 AFB90018 */  sw      $t9, 0x0018($sp)
/* 0005C 80B108CC AFB80014 */  sw      $t8, 0x0014($sp)
/* 00060 80B108D0 AFAF0010 */  sw      $t7, 0x0010($sp)
/* 00064 80B108D4 24E70338 */  addiu   $a3, $a3, 0x0338           ## $a3 = 06000338
/* 00068 80B108D8 24C69B38 */  addiu   $a2, $a2, 0x9B38           ## $a2 = 06009B38
/* 0006C 80B108DC 8FA40034 */  lw      $a0, 0x0034($sp)
/* 00070 80B108E0 0C0291BE */  jal     SkelAnime_InitSV
/* 00074 80B108E4 2605014C */  addiu   $a1, $s0, 0x014C           ## $a1 = 0000014C
/* 00078 80B108E8 3C088016 */  lui     $t0, 0x8016                ## $t0 = 80160000
/* 0007C 80B108EC 8D08E664 */  lw      $t0, -0x199C($t0)          ## 8015E664
/* 00080 80B108F0 3C0B80B1 */  lui     $t3, %hi(func_80B11310)    ## $t3 = 80B10000
/* 00084 80B108F4 240A0014 */  addiu   $t2, $zero, 0x0014         ## $t2 = 00000014
/* 00088 80B108F8 11000003 */  beq     $t0, $zero, .L80B10908
/* 0008C 80B108FC 256B1310 */  addiu   $t3, $t3, %lo(func_80B11310) ## $t3 = 80B11310
/* 00090 80B10900 24090014 */  addiu   $t1, $zero, 0x0014         ## $t1 = 00000014
/* 00094 80B10904 A6090204 */  sh      $t1, 0x0204($s0)           ## 00000204
.L80B10908:
/* 00098 80B10908 3C0D80B1 */  lui     $t5, %hi(func_80B10948)    ## $t5 = 80B10000
/* 0009C 80B1090C 240C0064 */  addiu   $t4, $zero, 0x0064         ## $t4 = 00000064
/* 000A0 80B10910 25AD0948 */  addiu   $t5, $t5, %lo(func_80B10948) ## $t5 = 80B10948
/* 000A4 80B10914 A60A020E */  sh      $t2, 0x020E($s0)           ## 0000020E
/* 000A8 80B10918 A600020C */  sh      $zero, 0x020C($s0)         ## 0000020C
/* 000AC 80B1091C AE0B0224 */  sw      $t3, 0x0224($s0)           ## 00000224
/* 000B0 80B10920 A60C00A8 */  sh      $t4, 0x00A8($s0)           ## 000000A8
/* 000B4 80B10924 AE0D01FC */  sw      $t5, 0x01FC($s0)           ## 000001FC
/* 000B8 80B10928 8FBF002C */  lw      $ra, 0x002C($sp)
/* 000BC 80B1092C 8FB00028 */  lw      $s0, 0x0028($sp)
/* 000C0 80B10930 27BD0030 */  addiu   $sp, $sp, 0x0030           ## $sp = 00000000
/* 000C4 80B10934 03E00008 */  jr      $ra
/* 000C8 80B10938 00000000 */  nop
