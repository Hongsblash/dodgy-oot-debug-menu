glabel func_80B02884
/* 00614 80B02884 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 00618 80B02888 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 0061C 80B0288C C4860080 */  lwc1    $f6, 0x0080($a0)           ## 00000080
/* 00620 80B02890 C4840540 */  lwc1    $f4, 0x0540($a0)           ## 00000540
/* 00624 80B02894 24050003 */  addiu   $a1, $zero, 0x0003         ## $a1 = 00000003
/* 00628 80B02898 46062200 */  add.s   $f8, $f4, $f6              
/* 0062C 80B0289C E4880028 */  swc1    $f8, 0x0028($a0)           ## 00000028
/* 00630 80B028A0 0C2C09C0 */  jal     func_80B02700              
/* 00634 80B028A4 AFA40018 */  sw      $a0, 0x0018($sp)           
/* 00638 80B028A8 4600028D */  trunc.w.s $f10, $f0                  
/* 0063C 80B028AC 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 00640 80B028B0 440F5000 */  mfc1    $t7, $f10                  
/* 00644 80B028B4 00000000 */  nop
/* 00648 80B028B8 A48F0534 */  sh      $t7, 0x0534($a0)           ## 00000534
/* 0064C 80B028BC 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 00650 80B028C0 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 00654 80B028C4 03E00008 */  jr      $ra                        
/* 00658 80B028C8 00000000 */  nop


