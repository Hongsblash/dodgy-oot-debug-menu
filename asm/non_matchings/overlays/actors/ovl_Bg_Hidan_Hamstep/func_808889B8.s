glabel func_808889B8
/* 00978 808889B8 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 0097C 808889BC AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 00980 808889C0 AFA50024 */  sw      $a1, 0x0024($sp)           
/* 00984 808889C4 8C820118 */  lw      $v0, 0x0118($a0)           ## 00000118
/* 00988 808889C8 AFA40020 */  sw      $a0, 0x0020($sp)           
/* 0098C 808889CC 0C2221CD */  jal     func_80888734              
/* 00990 808889D0 AFA20018 */  sw      $v0, 0x0018($sp)           
/* 00994 808889D4 8FA20018 */  lw      $v0, 0x0018($sp)           
/* 00998 808889D8 24010004 */  addiu   $at, $zero, 0x0004         ## $at = 00000004
/* 0099C 808889DC 8FA40020 */  lw      $a0, 0x0020($sp)           
/* 009A0 808889E0 8C430240 */  lw      $v1, 0x0240($v0)           ## 00000240
/* 009A4 808889E4 10610007 */  beq     $v1, $at, .L80888A04       
/* 009A8 808889E8 24010003 */  addiu   $at, $zero, 0x0003         ## $at = 00000003
/* 009AC 808889EC 54610017 */  bnel    $v1, $at, .L80888A4C       
/* 009B0 808889F0 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 009B4 808889F4 8C4E0244 */  lw      $t6, 0x0244($v0)           ## 00000244
/* 009B8 808889F8 29C10005 */  slti    $at, $t6, 0x0005           
/* 009BC 808889FC 54200013 */  bnel    $at, $zero, .L80888A4C     
/* 009C0 80888A00 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80888A04:
/* 009C4 80888A04 8482001C */  lh      $v0, 0x001C($a0)           ## 0000001C
/* 009C8 80888A08 24010001 */  addiu   $at, $zero, 0x0001         ## $at = 00000001
/* 009CC 80888A0C 304200FF */  andi    $v0, $v0, 0x00FF           ## $v0 = 00000000
/* 009D0 80888A10 1441000B */  bne     $v0, $at, .L80888A40       
/* 009D4 80888A14 00027880 */  sll     $t7, $v0,  2               
/* 009D8 80888A18 3C018089 */  lui     $at, %hi(.L80888D6C)       ## $at = 80890000
/* 009DC 80888A1C 002F0821 */  addu    $at, $at, $t7              
/* 009E0 80888A20 C4248D6C */  lwc1    $f4, %lo(.L80888D6C)($at)  
/* 009E4 80888A24 C486000C */  lwc1    $f6, 0x000C($a0)           ## 0000000C
/* 009E8 80888A28 24050004 */  addiu   $a1, $zero, 0x0004         ## $a1 = 00000004
/* 009EC 80888A2C 46062200 */  add.s   $f8, $f4, $f6              
/* 009F0 80888A30 0C222010 */  jal     func_80888040              
/* 009F4 80888A34 E4880028 */  swc1    $f8, 0x0028($a0)           ## 00000028
/* 009F8 80888A38 10000004 */  beq     $zero, $zero, .L80888A4C   
/* 009FC 80888A3C 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80888A40:
/* 00A00 80888A40 0C222010 */  jal     func_80888040              
/* 00A04 80888A44 24050003 */  addiu   $a1, $zero, 0x0003         ## $a1 = 00000003
/* 00A08 80888A48 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80888A4C:
/* 00A0C 80888A4C 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 00A10 80888A50 03E00008 */  jr      $ra                        
/* 00A14 80888A54 00000000 */  nop


