glabel func_80B50278
/* 018E8 80B50278 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 018EC 80B5027C 00803025 */  or      $a2, $a0, $zero            ## $a2 = 00000000
/* 018F0 80B50280 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 018F4 80B50284 00A02025 */  or      $a0, $a1, $zero            ## $a0 = 00000000
/* 018F8 80B50288 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 018FC 80B5028C 0C2D3B65 */  jal     func_80B4ED94              
/* 01900 80B50290 AFA60018 */  sw      $a2, 0x0018($sp)           
/* 01904 80B50294 8C4E000C */  lw      $t6, 0x000C($v0)           ## 0000000C
/* 01908 80B50298 8FA60018 */  lw      $a2, 0x0018($sp)           
/* 0190C 80B5029C 241900FF */  addiu   $t9, $zero, 0x00FF         ## $t9 = 000000FF
/* 01910 80B502A0 448E2000 */  mtc1    $t6, $f4                   ## $f4 = 0.00
/* 01914 80B502A4 24080002 */  addiu   $t0, $zero, 0x0002         ## $t0 = 00000002
/* 01918 80B502A8 24090001 */  addiu   $t1, $zero, 0x0001         ## $t1 = 00000001
/* 0191C 80B502AC 468021A0 */  cvt.s.w $f6, $f4                   
/* 01920 80B502B0 E4C60024 */  swc1    $f6, 0x0024($a2)           ## 00000024
/* 01924 80B502B4 8C4F0010 */  lw      $t7, 0x0010($v0)           ## 00000010
/* 01928 80B502B8 448F4000 */  mtc1    $t7, $f8                   ## $f8 = 0.00
/* 0192C 80B502BC 00000000 */  nop
/* 01930 80B502C0 468042A0 */  cvt.s.w $f10, $f8                  
/* 01934 80B502C4 E4CA0028 */  swc1    $f10, 0x0028($a2)          ## 00000028
/* 01938 80B502C8 8C580014 */  lw      $t8, 0x0014($v0)           ## 00000014
/* 0193C 80B502CC 44988000 */  mtc1    $t8, $f16                  ## $f16 = 0.00
/* 01940 80B502D0 00000000 */  nop
/* 01944 80B502D4 468084A0 */  cvt.s.w $f18, $f16                 
/* 01948 80B502D8 E4D2002C */  swc1    $f18, 0x002C($a2)          ## 0000002C
/* 0194C 80B502DC 84430008 */  lh      $v1, 0x0008($v0)           ## 00000008
/* 01950 80B502E0 A0D900C8 */  sb      $t9, 0x00C8($a2)           ## 000000C8
/* 01954 80B502E4 ACC8019C */  sw      $t0, 0x019C($a2)           ## 0000019C
/* 01958 80B502E8 ACC901A0 */  sw      $t1, 0x01A0($a2)           ## 000001A0
/* 0195C 80B502EC A4C300B6 */  sh      $v1, 0x00B6($a2)           ## 000000B6
/* 01960 80B502F0 A4C30032 */  sh      $v1, 0x0032($a2)           ## 00000032
/* 01964 80B502F4 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 01968 80B502F8 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 0196C 80B502FC 03E00008 */  jr      $ra                        
/* 01970 80B50300 00000000 */  nop


