glabel func_809C4E60
/* 00020 809C4E60 AFA50004 */  sw      $a1, 0x0004($sp)           
/* 00024 809C4E64 848E015C */  lh      $t6, 0x015C($a0)           ## 0000015C
/* 00028 809C4E68 3C18809C */  lui     $t8, %hi(func_809C4E8C)    ## $t8 = 809C0000
/* 0002C 809C4E6C 300F00FF */  andi    $t7, $zero, 0x00FF         ## $t7 = 00000000
/* 00030 809C4E70 11C00004 */  beq     $t6, $zero, .L809C4E84     
/* 00034 809C4E74 27184E8C */  addiu   $t8, $t8, %lo(func_809C4E8C) ## $t8 = 809C4E8C
/* 00038 809C4E78 A0800164 */  sb      $zero, 0x0164($a0)         ## 00000164
/* 0003C 809C4E7C A48F015C */  sh      $t7, 0x015C($a0)           ## 0000015C
/* 00040 809C4E80 AC98014C */  sw      $t8, 0x014C($a0)           ## 0000014C
.L809C4E84:
/* 00044 809C4E84 03E00008 */  jr      $ra                        
/* 00048 809C4E88 00000000 */  nop


