glabel func_80B060FC
/* 0186C 80B060FC C4800060 */  lwc1    $f0, 0x0060($a0)           ## 00000060
/* 01870 80B06100 C4840028 */  lwc1    $f4, 0x0028($a0)           ## 00000028
/* 01874 80B06104 C48A0080 */  lwc1    $f10, 0x0080($a0)          ## 00000080
/* 01878 80B06108 46000180 */  add.s   $f6, $f0, $f0              
/* 0187C 80B0610C C488040C */  lwc1    $f8, 0x040C($a0)           ## 0000040C
/* 01880 80B06110 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 01884 80B06114 46062080 */  add.s   $f2, $f4, $f6              
/* 01888 80B06118 460A1401 */  sub.s   $f16, $f2, $f10            
/* 0188C 80B0611C 4608803E */  c.le.s  $f16, $f8                  
/* 01890 80B06120 00000000 */  nop
/* 01894 80B06124 45000003 */  bc1f    .L80B06134                 
/* 01898 80B06128 00000000 */  nop
/* 0189C 80B0612C 03E00008 */  jr      $ra                        
/* 018A0 80B06130 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
.L80B06134:
/* 018A4 80B06134 03E00008 */  jr      $ra                        
/* 018A8 80B06138 00000000 */  nop


