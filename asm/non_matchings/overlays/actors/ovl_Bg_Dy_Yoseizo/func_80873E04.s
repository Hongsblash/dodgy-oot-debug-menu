glabel func_80873E04
/* 015D4 80873E04 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 015D8 80873E08 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 015DC 80873E0C AFA40018 */  sw      $a0, 0x0018($sp)           
/* 015E0 80873E10 848F02E8 */  lh      $t7, 0x02E8($a0)           ## 000002E8
/* 015E4 80873E14 00807025 */  or      $t6, $a0, $zero            ## $t6 = 00000000
/* 015E8 80873E18 00A03825 */  or      $a3, $a1, $zero            ## $a3 = 00000000
/* 015EC 80873E1C 15E0001D */  bne     $t7, $zero, .L80873E94     
/* 015F0 80873E20 00A02025 */  or      $a0, $a1, $zero            ## $a0 = 00000000
/* 015F4 80873E24 01C02825 */  or      $a1, $t6, $zero            ## $a1 = 00000000
/* 015F8 80873E28 24060007 */  addiu   $a2, $zero, 0x0007         ## $a2 = 00000007
/* 015FC 80873E2C 0C00B7D5 */  jal     func_8002DF54              
/* 01600 80873E30 AFA7001C */  sw      $a3, 0x001C($sp)           
/* 01604 80873E34 8FA7001C */  lw      $a3, 0x001C($sp)           
/* 01608 80873E38 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 0160C 80873E3C 2402012E */  addiu   $v0, $zero, 0x012E         ## $v0 = 0000012E
/* 01610 80873E40 00270821 */  addu    $at, $at, $a3              
/* 01614 80873E44 A0200AE3 */  sb      $zero, 0x0AE3($at)         ## 00010AE3
/* 01618 80873E48 8CE41C64 */  lw      $a0, 0x1C64($a3)           ## 00001C64
/* 0161C 80873E4C 5080000D */  beql    $a0, $zero, .L80873E84     
/* 01620 80873E50 00E02025 */  or      $a0, $a3, $zero            ## $a0 = 00000000
/* 01624 80873E54 84980000 */  lh      $t8, 0x0000($a0)           ## 00000000
.L80873E58:
/* 01628 80873E58 10580003 */  beq     $v0, $t8, .L80873E68       
/* 0162C 80873E5C 00000000 */  nop
/* 01630 80873E60 10000005 */  beq     $zero, $zero, .L80873E78   
/* 01634 80873E64 8C840124 */  lw      $a0, 0x0124($a0)           ## 00000124
.L80873E68:
/* 01638 80873E68 0C00B55C */  jal     Actor_Kill
              
/* 0163C 80873E6C AFA7001C */  sw      $a3, 0x001C($sp)           
/* 01640 80873E70 10000003 */  beq     $zero, $zero, .L80873E80   
/* 01644 80873E74 8FA7001C */  lw      $a3, 0x001C($sp)           
.L80873E78:
/* 01648 80873E78 5480FFF7 */  bnel    $a0, $zero, .L80873E58     
/* 0164C 80873E7C 84980000 */  lh      $t8, 0x0000($a0)           ## 00000000
.L80873E80:
/* 01650 80873E80 00E02025 */  or      $a0, $a3, $zero            ## $a0 = 00000000
.L80873E84:
/* 01654 80873E84 0C00B2ED */  jal     Flags_UnsetSwitch
              
/* 01658 80873E88 24050038 */  addiu   $a1, $zero, 0x0038         ## $a1 = 00000038
/* 0165C 80873E8C 0C00B55C */  jal     Actor_Kill
              
/* 01660 80873E90 8FA40018 */  lw      $a0, 0x0018($sp)           
.L80873E94:
/* 01664 80873E94 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 01668 80873E98 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 0166C 80873E9C 03E00008 */  jr      $ra                        
/* 01670 80873EA0 00000000 */  nop


