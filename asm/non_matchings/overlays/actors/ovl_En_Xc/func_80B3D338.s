glabel func_80B3D338
/* 01158 80B3D338 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 0115C 80B3D33C AFA40018 */  sw      $a0, 0x0018($sp)           
/* 01160 80B3D340 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 01164 80B3D344 00A02025 */  or      $a0, $a1, $zero            ## $a0 = 00000000
/* 01168 80B3D348 0C2CF134 */  jal     func_80B3C4D0              
/* 0116C 80B3D34C 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 01170 80B3D350 8FA50018 */  lw      $a1, 0x0018($sp)           
/* 01174 80B3D354 8CA40310 */  lw      $a0, 0x0310($a1)           ## 00000310
/* 01178 80B3D358 50800012 */  beql    $a0, $zero, .L80B3D3A4     
/* 0117C 80B3D35C 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 01180 80B3D360 8C4E000C */  lw      $t6, 0x000C($v0)           ## 0000000C
/* 01184 80B3D364 24830024 */  addiu   $v1, $a0, 0x0024           ## $v1 = 00000024
/* 01188 80B3D368 448E2000 */  mtc1    $t6, $f4                   ## $f4 = 0.00
/* 0118C 80B3D36C 00000000 */  nop
/* 01190 80B3D370 468021A0 */  cvt.s.w $f6, $f4                   
/* 01194 80B3D374 E4660000 */  swc1    $f6, 0x0000($v1)           ## 00000024
/* 01198 80B3D378 8C4F0010 */  lw      $t7, 0x0010($v0)           ## 00000010
/* 0119C 80B3D37C 448F4000 */  mtc1    $t7, $f8                   ## $f8 = 0.00
/* 011A0 80B3D380 00000000 */  nop
/* 011A4 80B3D384 468042A0 */  cvt.s.w $f10, $f8                  
/* 011A8 80B3D388 E46A0004 */  swc1    $f10, 0x0004($v1)          ## 00000028
/* 011AC 80B3D38C 8C580014 */  lw      $t8, 0x0014($v0)           ## 00000014
/* 011B0 80B3D390 44988000 */  mtc1    $t8, $f16                  ## $f16 = 0.00
/* 011B4 80B3D394 00000000 */  nop
/* 011B8 80B3D398 468084A0 */  cvt.s.w $f18, $f16                 
/* 011BC 80B3D39C E4720008 */  swc1    $f18, 0x0008($v1)          ## 0000002C
/* 011C0 80B3D3A0 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80B3D3A4:
/* 011C4 80B3D3A4 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 011C8 80B3D3A8 03E00008 */  jr      $ra                        
/* 011CC 80B3D3AC 00000000 */  nop


