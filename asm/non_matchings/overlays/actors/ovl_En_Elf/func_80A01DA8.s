glabel func_80A01DA8
/* 00178 80A01DA8 44801000 */  mtc1    $zero, $f2                 ## $f2 = 0.00
/* 0017C 80A01DAC 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 00180 80A01DB0 44810000 */  mtc1    $at, $f0                   ## $f0 = 1.00
/* 00184 80A01DB4 3C0D80A0 */  lui     $t5, %hi(func_80A02A20)    ## $t5 = 80A00000
/* 00188 80A01DB8 240B1000 */  addiu   $t3, $zero, 0x1000         ## $t3 = 00001000
/* 0018C 80A01DBC 240C0200 */  addiu   $t4, $zero, 0x0200         ## $t4 = 00000200
/* 00190 80A01DC0 25AD2A20 */  addiu   $t5, $t5, %lo(func_80A02A20) ## $t5 = 80A02A20
/* 00194 80A01DC4 A48B02AE */  sh      $t3, 0x02AE($a0)           ## 000002AE
/* 00198 80A01DC8 A48C02B0 */  sh      $t4, 0x02B0($a0)           ## 000002B0
/* 0019C 80A01DCC AC8D02C8 */  sw      $t5, 0x02C8($a0)           ## 000002C8
/* 001A0 80A01DD0 E48202B4 */  swc1    $f2, 0x02B4($a0)           ## 000002B4
/* 001A4 80A01DD4 E48202B8 */  swc1    $f2, 0x02B8($a0)           ## 000002B8
/* 001A8 80A01DD8 03E00008 */  jr      $ra                        
/* 001AC 80A01DDC E4800168 */  swc1    $f0, 0x0168($a0)           ## 00000168


