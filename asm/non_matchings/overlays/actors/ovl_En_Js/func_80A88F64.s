glabel func_80A88F64
/* 00154 80A88F64 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 00158 80A88F68 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 0015C 80A88F6C AFA5001C */  sw      $a1, 0x001C($sp)           
/* 00160 80A88F70 AFA60020 */  sw      $a2, 0x0020($sp)           
/* 00164 80A88F74 0C00BC65 */  jal     func_8002F194              
/* 00168 80A88F78 AFA40018 */  sw      $a0, 0x0018($sp)           
/* 0016C 80A88F7C 10400003 */  beq     $v0, $zero, .L80A88F8C     
/* 00170 80A88F80 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 00174 80A88F84 1000001C */  beq     $zero, $zero, .L80A88FF8   
/* 00178 80A88F88 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
.L80A88F8C:
/* 0017C 80A88F8C 848F008A */  lh      $t7, 0x008A($a0)           ## 0000008A
/* 00180 80A88F90 849800B6 */  lh      $t8, 0x00B6($a0)           ## 000000B6
/* 00184 80A88F94 97AE0022 */  lhu     $t6, 0x0022($sp)           
/* 00188 80A88F98 01F81023 */  subu    $v0, $t7, $t8              
/* 0018C 80A88F9C 00021400 */  sll     $v0, $v0, 16               
/* 00190 80A88FA0 00021403 */  sra     $v0, $v0, 16               
/* 00194 80A88FA4 04400003 */  bltz    $v0, .L80A88FB4            
/* 00198 80A88FA8 A48E010E */  sh      $t6, 0x010E($a0)           ## 0000010E
/* 0019C 80A88FAC 10000002 */  beq     $zero, $zero, .L80A88FB8   
/* 001A0 80A88FB0 00401825 */  or      $v1, $v0, $zero            ## $v1 = 00000001
.L80A88FB4:
/* 001A4 80A88FB4 00021823 */  subu    $v1, $zero, $v0            
.L80A88FB8:
/* 001A8 80A88FB8 28611801 */  slti    $at, $v1, 0x1801           
/* 001AC 80A88FBC 1020000D */  beq     $at, $zero, .L80A88FF4     
/* 001B0 80A88FC0 3C0142C8 */  lui     $at, 0x42C8                ## $at = 42C80000
/* 001B4 80A88FC4 44810000 */  mtc1    $at, $f0                   ## $f0 = 100.00
/* 001B8 80A88FC8 C4840090 */  lwc1    $f4, 0x0090($a0)           ## 00000090
/* 001BC 80A88FCC 4600203C */  c.lt.s  $f4, $f0                   
/* 001C0 80A88FD0 00000000 */  nop
/* 001C4 80A88FD4 45020008 */  bc1fl   .L80A88FF8                 
/* 001C8 80A88FD8 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 001CC 80A88FDC 94990284 */  lhu     $t9, 0x0284($a0)           ## 00000284
/* 001D0 80A88FE0 44060000 */  mfc1    $a2, $f0                   
/* 001D4 80A88FE4 37280001 */  ori     $t0, $t9, 0x0001           ## $t0 = 00000001
/* 001D8 80A88FE8 A4880284 */  sh      $t0, 0x0284($a0)           ## 00000284
/* 001DC 80A88FEC 0C00BCB3 */  jal     func_8002F2CC              
/* 001E0 80A88FF0 8FA5001C */  lw      $a1, 0x001C($sp)           
.L80A88FF4:
/* 001E4 80A88FF4 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
.L80A88FF8:
/* 001E8 80A88FF8 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 001EC 80A88FFC 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 001F0 80A89000 03E00008 */  jr      $ra                        
/* 001F4 80A89004 00000000 */  nop
