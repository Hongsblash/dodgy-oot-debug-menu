glabel func_80954924
/* 05664 80954924 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 05668 80954928 AFA40028 */  sw      $a0, 0x0028($sp)
/* 0566C 8095492C AFBF0024 */  sw      $ra, 0x0024($sp)
/* 05670 80954930 3C040602 */  lui     $a0, 0x0602                ## $a0 = 06020000
/* 05674 80954934 AFA5002C */  sw      $a1, 0x002C($sp)
/* 05678 80954938 0C028800 */  jal     SkelAnime_GetFrameCount

/* 0567C 8095493C 24848150 */  addiu   $a0, $a0, 0x8150           ## $a0 = 06018150
/* 05680 80954940 44822000 */  mtc1    $v0, $f4                   ## $f4 = 0.00
/* 05684 80954944 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 05688 80954948 8FA40028 */  lw      $a0, 0x0028($sp)
/* 0568C 8095494C 468021A0 */  cvt.s.w $f6, $f4
/* 05690 80954950 3C050602 */  lui     $a1, 0x0602                ## $a1 = 06020000
/* 05694 80954954 240E0002 */  addiu   $t6, $zero, 0x0002         ## $t6 = 00000002
/* 05698 80954958 44070000 */  mfc1    $a3, $f0
/* 0569C 8095495C AFAE0014 */  sw      $t6, 0x0014($sp)
/* 056A0 80954960 24A58150 */  addiu   $a1, $a1, 0x8150           ## $a1 = 06018150
/* 056A4 80954964 E7A60010 */  swc1    $f6, 0x0010($sp)
/* 056A8 80954968 3C063F80 */  lui     $a2, 0x3F80                ## $a2 = 3F800000
/* 056AC 8095496C 2484014C */  addiu   $a0, $a0, 0x014C           ## $a0 = 0000014C
/* 056B0 80954970 0C029468 */  jal     SkelAnime_ChangeAnim

/* 056B4 80954974 E7A00018 */  swc1    $f0, 0x0018($sp)
/* 056B8 80954978 8FA40028 */  lw      $a0, 0x0028($sp)
/* 056BC 8095497C 2401FFFE */  addiu   $at, $zero, 0xFFFE         ## $at = FFFFFFFE
/* 056C0 80954980 3C058095 */  lui     $a1, %hi(func_809549A8)    ## $a1 = 80950000
/* 056C4 80954984 8C8F0004 */  lw      $t7, 0x0004($a0)           ## 00000004
/* 056C8 80954988 24A549A8 */  addiu   $a1, $a1, %lo(func_809549A8) ## $a1 = 809549A8
/* 056CC 8095498C 01E1C024 */  and     $t8, $t7, $at
/* 056D0 80954990 0C253CB0 */  jal     func_8094F2C0
/* 056D4 80954994 AC980004 */  sw      $t8, 0x0004($a0)           ## 00000004
/* 056D8 80954998 8FBF0024 */  lw      $ra, 0x0024($sp)
/* 056DC 8095499C 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 056E0 809549A0 03E00008 */  jr      $ra
/* 056E4 809549A4 00000000 */  nop
