glabel func_80A59828
/* 00C18 80A59828 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 00C1C 80A5982C AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 00C20 80A59830 AFA40020 */  sw      $a0, 0x0020($sp)           
/* 00C24 80A59834 0C296306 */  jal     func_80A58C18              
/* 00C28 80A59838 AFA50024 */  sw      $a1, 0x0024($sp)           
/* 00C2C 80A5983C 14400013 */  bne     $v0, $zero, .L80A5988C     
/* 00C30 80A59840 8FA60024 */  lw      $a2, 0x0024($sp)           
/* 00C34 80A59844 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 00C38 80A59848 00C11021 */  addu    $v0, $a2, $at              
/* 00C3C 80A5984C 804E1CED */  lb      $t6, 0x1CED($v0)           ## 00001CED
/* 00C40 80A59850 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 00C44 80A59854 34211CBC */  ori     $at, $at, 0x1CBC           ## $at = 00011CBC
/* 00C48 80A59858 15C0000C */  bne     $t6, $zero, .L80A5988C     
/* 00C4C 80A5985C 00C02025 */  or      $a0, $a2, $zero            ## $a0 = 00000000
/* 00C50 80A59860 00C12821 */  addu    $a1, $a2, $at              
/* 00C54 80A59864 0C025D4D */  jal     func_80097534              
/* 00C58 80A59868 AFA2001C */  sw      $v0, 0x001C($sp)           
/* 00C5C 80A5986C 8FA2001C */  lw      $v0, 0x001C($sp)           
/* 00C60 80A59870 8FB80020 */  lw      $t8, 0x0020($sp)           
/* 00C64 80A59874 844F1E18 */  lh      $t7, 0x1E18($v0)           ## 00001E18
/* 00C68 80A59878 15E00002 */  bne     $t7, $zero, .L80A59884     
/* 00C6C 80A5987C 00000000 */  nop
/* 00C70 80A59880 A300014F */  sb      $zero, 0x014F($t8)         ## 0000014F
.L80A59884:
/* 00C74 80A59884 0C296312 */  jal     func_80A58C48              
/* 00C78 80A59888 8FA40020 */  lw      $a0, 0x0020($sp)           
.L80A5988C:
/* 00C7C 80A5988C 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 00C80 80A59890 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 00C84 80A59894 03E00008 */  jr      $ra                        
/* 00C88 80A59898 00000000 */  nop


