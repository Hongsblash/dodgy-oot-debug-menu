glabel func_80A610C8
/* 05DD8 80A610C8 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 05DDC 80A610CC AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 05DE0 80A610D0 90AE1D6C */  lbu     $t6, 0x1D6C($a1)           ## 00001D6C
/* 05DE4 80A610D4 24010003 */  addiu   $at, $zero, 0x0003         ## $at = 00000003
/* 05DE8 80A610D8 00803825 */  or      $a3, $a0, $zero            ## $a3 = 00000000
/* 05DEC 80A610DC 15C1000A */  bne     $t6, $at, .L80A61108       
/* 05DF0 80A610E0 8CA61D88 */  lw      $a2, 0x1D88($a1)           ## 00001D88
/* 05DF4 80A610E4 240F0001 */  addiu   $t7, $zero, 0x0001         ## $t7 = 00000001
/* 05DF8 80A610E8 2418000A */  addiu   $t8, $zero, 0x000A         ## $t8 = 0000000A
/* 05DFC 80A610EC 24190002 */  addiu   $t9, $zero, 0x0002         ## $t9 = 00000002
/* 05E00 80A610F0 AC8F020C */  sw      $t7, 0x020C($a0)           ## 0000020C
/* 05E04 80A610F4 A498001C */  sh      $t8, 0x001C($a0)           ## 0000001C
/* 05E08 80A610F8 0C29723F */  jal     func_80A5C8FC              
/* 05E0C 80A610FC AC99014C */  sw      $t9, 0x014C($a0)           ## 0000014C
/* 05E10 80A61100 10000043 */  beq     $zero, $zero, .L80A61210   
/* 05E14 80A61104 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80A61108:
/* 05E18 80A61108 50C00041 */  beql    $a2, $zero, .L80A61210     
/* 05E1C 80A6110C 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 05E20 80A61110 94C40000 */  lhu     $a0, 0x0000($a2)           ## 00000000
/* 05E24 80A61114 AFA70020 */  sw      $a3, 0x0020($sp)           
/* 05E28 80A61118 AFA60018 */  sw      $a2, 0x0018($sp)           
/* 05E2C 80A6111C 0C29841E */  jal     func_80A61078              
/* 05E30 80A61120 AFA50024 */  sw      $a1, 0x0024($sp)           
/* 05E34 80A61124 8FA50024 */  lw      $a1, 0x0024($sp)           
/* 05E38 80A61128 8FA60018 */  lw      $a2, 0x0018($sp)           
/* 05E3C 80A6112C 8FA70020 */  lw      $a3, 0x0020($sp)           
/* 05E40 80A61130 10400036 */  beq     $v0, $zero, .L80A6120C     
/* 05E44 80A61134 00404025 */  or      $t0, $v0, $zero            ## $t0 = 00000000
/* 05E48 80A61138 8CE30380 */  lw      $v1, 0x0380($a3)           ## 00000380
/* 05E4C 80A6113C 00E02025 */  or      $a0, $a3, $zero            ## $a0 = 00000000
/* 05E50 80A61140 0008C880 */  sll     $t9, $t0,  2               
/* 05E54 80A61144 1043002A */  beq     $v0, $v1, .L80A611F0       
/* 05E58 80A61148 3C0180A6 */  lui     $at, %hi(D_80A666A4)       ## $at = 80A60000
/* 05E5C 80A6114C 1460001E */  bne     $v1, $zero, .L80A611C8     
/* 05E60 80A61150 00390821 */  addu    $at, $at, $t9              
/* 05E64 80A61154 8CC9000C */  lw      $t1, 0x000C($a2)           ## 0000000C
/* 05E68 80A61158 44892000 */  mtc1    $t1, $f4                   ## $f4 = 0.00
/* 05E6C 80A6115C 00000000 */  nop
/* 05E70 80A61160 468021A0 */  cvt.s.w $f6, $f4                   
/* 05E74 80A61164 E4E60024 */  swc1    $f6, 0x0024($a3)           ## 00000024
/* 05E78 80A61168 8CCA0010 */  lw      $t2, 0x0010($a2)           ## 00000010
/* 05E7C 80A6116C 8CF80024 */  lw      $t8, 0x0024($a3)           ## 00000024
/* 05E80 80A61170 448A4000 */  mtc1    $t2, $f8                   ## $f8 = 0.00
/* 05E84 80A61174 00000000 */  nop
/* 05E88 80A61178 468042A0 */  cvt.s.w $f10, $f8                  
/* 05E8C 80A6117C E4EA0028 */  swc1    $f10, 0x0028($a3)          ## 00000028
/* 05E90 80A61180 8CCB0014 */  lw      $t3, 0x0014($a2)           ## 00000014
/* 05E94 80A61184 8CEF0028 */  lw      $t7, 0x0028($a3)           ## 00000028
/* 05E98 80A61188 448B8000 */  mtc1    $t3, $f16                  ## $f16 = 0.00
/* 05E9C 80A6118C 00000000 */  nop
/* 05EA0 80A61190 468084A0 */  cvt.s.w $f18, $f16                 
/* 05EA4 80A61194 E4F2002C */  swc1    $f18, 0x002C($a3)          ## 0000002C
/* 05EA8 80A61198 94CC0008 */  lhu     $t4, 0x0008($a2)           ## 00000008
/* 05EAC 80A6119C ACF80100 */  sw      $t8, 0x0100($a3)           ## 00000100
/* 05EB0 80A611A0 8CF8002C */  lw      $t8, 0x002C($a3)           ## 0000002C
/* 05EB4 80A611A4 A4EC0032 */  sh      $t4, 0x0032($a3)           ## 00000032
/* 05EB8 80A611A8 88EE0030 */  lwl     $t6, 0x0030($a3)           ## 00000030
/* 05EBC 80A611AC 98EE0033 */  lwr     $t6, 0x0033($a3)           ## 00000033
/* 05EC0 80A611B0 ACEF0104 */  sw      $t7, 0x0104($a3)           ## 00000104
/* 05EC4 80A611B4 ACF80108 */  sw      $t8, 0x0108($a3)           ## 00000108
/* 05EC8 80A611B8 A8EE00B4 */  swl     $t6, 0x00B4($a3)           ## 000000B4
/* 05ECC 80A611BC B8EE00B7 */  swr     $t6, 0x00B7($a3)           ## 000000B7
/* 05ED0 80A611C0 94EE0034 */  lhu     $t6, 0x0034($a3)           ## 00000034
/* 05ED4 80A611C4 A4EE00B8 */  sh      $t6, 0x00B8($a3)           ## 000000B8
.L80A611C8:
/* 05ED8 80A611C8 ACE80380 */  sw      $t0, 0x0380($a3)           ## 00000380
/* 05EDC 80A611CC 8C3966A4 */  lw      $t9, %lo(D_80A666A4)($at)  
/* 05EE0 80A611D0 AFA70020 */  sw      $a3, 0x0020($sp)           
/* 05EE4 80A611D4 AFA60018 */  sw      $a2, 0x0018($sp)           
/* 05EE8 80A611D8 0320F809 */  jalr    $ra, $t9                   
/* 05EEC 80A611DC AFA50024 */  sw      $a1, 0x0024($sp)           
/* 05EF0 80A611E0 8FA70020 */  lw      $a3, 0x0020($sp)           
/* 05EF4 80A611E4 8FA50024 */  lw      $a1, 0x0024($sp)           
/* 05EF8 80A611E8 8FA60018 */  lw      $a2, 0x0018($sp)           
/* 05EFC 80A611EC 8CE30380 */  lw      $v1, 0x0380($a3)           ## 00000380
.L80A611F0:
/* 05F00 80A611F0 00034880 */  sll     $t1, $v1,  2               
/* 05F04 80A611F4 3C1980A6 */  lui     $t9, %hi(D_80A666BC)       ## $t9 = 80A60000
/* 05F08 80A611F8 0329C821 */  addu    $t9, $t9, $t1              
/* 05F0C 80A611FC 8F3966BC */  lw      $t9, %lo(D_80A666BC)($t9)  
/* 05F10 80A61200 00E02025 */  or      $a0, $a3, $zero            ## $a0 = 00000000
/* 05F14 80A61204 0320F809 */  jalr    $ra, $t9                   
/* 05F18 80A61208 00000000 */  nop
.L80A6120C:
/* 05F1C 80A6120C 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80A61210:
/* 05F20 80A61210 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 05F24 80A61214 03E00008 */  jr      $ra                        
/* 05F28 80A61218 00000000 */  nop


