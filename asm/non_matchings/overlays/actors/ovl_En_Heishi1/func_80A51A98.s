glabel func_80A51A98
/* 007C8 80A51A98 27BDFFD0 */  addiu   $sp, $sp, 0xFFD0           ## $sp = FFFFFFD0
/* 007CC 80A51A9C AFB00028 */  sw      $s0, 0x0028($sp)           
/* 007D0 80A51AA0 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 007D4 80A51AA4 AFBF002C */  sw      $ra, 0x002C($sp)           
/* 007D8 80A51AA8 3C040600 */  lui     $a0, 0x0600                ## $a0 = 06000000
/* 007DC 80A51AAC AFA50034 */  sw      $a1, 0x0034($sp)           
/* 007E0 80A51AB0 0C028800 */  jal     SkelAnime_GetFrameCount
              
/* 007E4 80A51AB4 24845880 */  addiu   $a0, $a0, 0x5880           ## $a0 = 06005880
/* 007E8 80A51AB8 44822000 */  mtc1    $v0, $f4                   ## $f4 = 0.00
/* 007EC 80A51ABC 3C01C040 */  lui     $at, 0xC040                ## $at = C0400000
/* 007F0 80A51AC0 44819000 */  mtc1    $at, $f18                  ## $f18 = -3.00
/* 007F4 80A51AC4 468021A0 */  cvt.s.w $f6, $f4                   
/* 007F8 80A51AC8 3C050600 */  lui     $a1, 0x0600                ## $a1 = 06000000
/* 007FC 80A51ACC 24A55880 */  addiu   $a1, $a1, 0x5880           ## $a1 = 06005880
/* 00800 80A51AD0 2604014C */  addiu   $a0, $s0, 0x014C           ## $a0 = 0000014C
/* 00804 80A51AD4 3C064040 */  lui     $a2, 0x4040                ## $a2 = 40400000
/* 00808 80A51AD8 24070000 */  addiu   $a3, $zero, 0x0000         ## $a3 = 00000000
/* 0080C 80A51ADC 4600320D */  trunc.w.s $f8, $f6                   
/* 00810 80A51AE0 AFA00014 */  sw      $zero, 0x0014($sp)         
/* 00814 80A51AE4 E7B20018 */  swc1    $f18, 0x0018($sp)          
/* 00818 80A51AE8 440F4000 */  mfc1    $t7, $f8                   
/* 0081C 80A51AEC 00000000 */  nop
/* 00820 80A51AF0 000FC400 */  sll     $t8, $t7, 16               
/* 00824 80A51AF4 0018CC03 */  sra     $t9, $t8, 16               
/* 00828 80A51AF8 44995000 */  mtc1    $t9, $f10                  ## $f10 = 0.00
/* 0082C 80A51AFC 00000000 */  nop
/* 00830 80A51B00 46805420 */  cvt.s.w $f16, $f10                 
/* 00834 80A51B04 0C029468 */  jal     SkelAnime_ChangeAnimation
              
/* 00838 80A51B08 E7B00010 */  swc1    $f16, 0x0010($sp)          
/* 0083C 80A51B0C 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 00840 80A51B10 2405702D */  addiu   $a1, $zero, 0x702D         ## $a1 = 0000702D
/* 00844 80A51B14 02003025 */  or      $a2, $s0, $zero            ## $a2 = 00000000
/* 00848 80A51B18 E6000278 */  swc1    $f0, 0x0278($s0)           ## 00000278
/* 0084C 80A51B1C E6000274 */  swc1    $f0, 0x0274($s0)           ## 00000274
/* 00850 80A51B20 0C042DA0 */  jal     func_8010B680              
/* 00854 80A51B24 8FA40034 */  lw      $a0, 0x0034($sp)           
/* 00858 80A51B28 8FA40034 */  lw      $a0, 0x0034($sp)           
/* 0085C 80A51B2C 0C021BC0 */  jal     Interface_SetDoAction              
/* 00860 80A51B30 24050012 */  addiu   $a1, $zero, 0x0012         ## $a1 = 00000012
/* 00864 80A51B34 3C0880A5 */  lui     $t0, %hi(func_80A51B54)    ## $t0 = 80A50000
/* 00868 80A51B38 25081B54 */  addiu   $t0, $t0, %lo(func_80A51B54) ## $t0 = 80A51B54
/* 0086C 80A51B3C AE08025C */  sw      $t0, 0x025C($s0)           ## 0000025C
/* 00870 80A51B40 8FBF002C */  lw      $ra, 0x002C($sp)           
/* 00874 80A51B44 8FB00028 */  lw      $s0, 0x0028($sp)           
/* 00878 80A51B48 27BD0030 */  addiu   $sp, $sp, 0x0030           ## $sp = 00000000
/* 0087C 80A51B4C 03E00008 */  jr      $ra                        
/* 00880 80A51B50 00000000 */  nop


