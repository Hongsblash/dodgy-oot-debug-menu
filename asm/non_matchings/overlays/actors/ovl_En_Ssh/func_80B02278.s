glabel func_80B02278
/* 00008 80B02278 27BDFFC0 */  addiu   $sp, $sp, 0xFFC0           ## $sp = FFFFFFC0
/* 0000C 80B0227C 3C0E80B0 */  lui     $t6, %hi(D_80B04508)       ## $t6 = 80B00000
/* 00010 80B02280 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 00014 80B02284 AFA40040 */  sw      $a0, 0x0040($sp)           
/* 00018 80B02288 AFA50044 */  sw      $a1, 0x0044($sp)           
/* 0001C 80B0228C 25CE4508 */  addiu   $t6, $t6, %lo(D_80B04508)  ## $t6 = 80B04508
/* 00020 80B02290 8DD80000 */  lw      $t8, 0x0000($t6)           ## 80B04508
/* 00024 80B02294 27A60034 */  addiu   $a2, $sp, 0x0034           ## $a2 = FFFFFFF4
/* 00028 80B02298 24080064 */  addiu   $t0, $zero, 0x0064         ## $t0 = 00000064
/* 0002C 80B0229C ACD80000 */  sw      $t8, 0x0000($a2)           ## FFFFFFF4
/* 00030 80B022A0 8DCF0004 */  lw      $t7, 0x0004($t6)           ## 80B0450C
/* 00034 80B022A4 240900DC */  addiu   $t1, $zero, 0x00DC         ## $t1 = 000000DC
/* 00038 80B022A8 240A0008 */  addiu   $t2, $zero, 0x0008         ## $t2 = 00000008
/* 0003C 80B022AC ACCF0004 */  sw      $t7, 0x0004($a2)           ## FFFFFFF8
/* 00040 80B022B0 8DD80008 */  lw      $t8, 0x0008($t6)           ## 80B04510
/* 00044 80B022B4 00C03825 */  or      $a3, $a2, $zero            ## $a3 = FFFFFFF4
/* 00048 80B022B8 27A50028 */  addiu   $a1, $sp, 0x0028           ## $a1 = FFFFFFE8
/* 0004C 80B022BC ACD80008 */  sw      $t8, 0x0008($a2)           ## FFFFFFFC
/* 00050 80B022C0 8FB90040 */  lw      $t9, 0x0040($sp)           
/* 00054 80B022C4 8FA40044 */  lw      $a0, 0x0044($sp)           
/* 00058 80B022C8 C7240024 */  lwc1    $f4, 0x0024($t9)           ## 00000024
/* 0005C 80B022CC E7A40028 */  swc1    $f4, 0x0028($sp)           
/* 00060 80B022D0 C7260080 */  lwc1    $f6, 0x0080($t9)           ## 00000080
/* 00064 80B022D4 E7A6002C */  swc1    $f6, 0x002C($sp)           
/* 00068 80B022D8 C728002C */  lwc1    $f8, 0x002C($t9)           ## 0000002C
/* 0006C 80B022DC AFAA0018 */  sw      $t2, 0x0018($sp)           
/* 00070 80B022E0 AFA90014 */  sw      $t1, 0x0014($sp)           
/* 00074 80B022E4 AFA80010 */  sw      $t0, 0x0010($sp)           
/* 00078 80B022E8 0C00A3E1 */  jal     func_80028F84              
/* 0007C 80B022EC E7A80030 */  swc1    $f8, 0x0030($sp)           
/* 00080 80B022F0 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 00084 80B022F4 27BD0040 */  addiu   $sp, $sp, 0x0040           ## $sp = 00000000
/* 00088 80B022F8 03E00008 */  jr      $ra                        
/* 0008C 80B022FC 00000000 */  nop


