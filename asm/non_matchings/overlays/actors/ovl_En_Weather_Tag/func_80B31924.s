glabel func_80B31924
/* 00904 80B31924 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 00908 80B31928 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 0090C 80B3192C 240E0002 */  addiu   $t6, $zero, 0x0002         ## $t6 = 00000002
/* 00910 80B31930 240F0064 */  addiu   $t7, $zero, 0x0064         ## $t7 = 00000064
/* 00914 80B31934 24180004 */  addiu   $t8, $zero, 0x0004         ## $t8 = 00000004
/* 00918 80B31938 AFA40028 */  sw      $a0, 0x0028($sp)           
/* 0091C 80B3193C AFA5002C */  sw      $a1, 0x002C($sp)           
/* 00920 80B31940 AFB8001C */  sw      $t8, 0x001C($sp)           
/* 00924 80B31944 AFAF0018 */  sw      $t7, 0x0018($sp)           
/* 00928 80B31948 AFAE0014 */  sw      $t6, 0x0014($sp)           
/* 0092C 80B3194C AFA00010 */  sw      $zero, 0x0010($sp)         
/* 00930 80B31950 00003025 */  or      $a2, $zero, $zero          ## $a2 = 00000000
/* 00934 80B31954 0C2CC4B2 */  jal     func_80B312C8              
/* 00938 80B31958 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 0093C 80B3195C 5040000D */  beql    $v0, $zero, .L80B31994     
/* 00940 80B31960 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 00944 80B31964 0C01DD89 */  jal     func_80077624              
/* 00948 80B31968 8FA4002C */  lw      $a0, 0x002C($sp)           
/* 0094C 80B3196C 8FA8002C */  lw      $t0, 0x002C($sp)           
/* 00950 80B31970 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 00954 80B31974 24190019 */  addiu   $t9, $zero, 0x0019         ## $t9 = 00000019
/* 00958 80B31978 00280821 */  addu    $at, $at, $t0              
/* 0095C 80B3197C A0390B12 */  sb      $t9, 0x0B12($at)           ## 00010B12
/* 00960 80B31980 3C0580B3 */  lui     $a1, %hi(func_80B319A0)    ## $a1 = 80B30000
/* 00964 80B31984 24A519A0 */  addiu   $a1, $a1, %lo(func_80B319A0) ## $a1 = 80B319A0
/* 00968 80B31988 0C2CC408 */  jal     func_80B31020              
/* 0096C 80B3198C 8FA40028 */  lw      $a0, 0x0028($sp)           
/* 00970 80B31990 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80B31994:
/* 00974 80B31994 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 00978 80B31998 03E00008 */  jr      $ra                        
/* 0097C 80B3199C 00000000 */  nop
