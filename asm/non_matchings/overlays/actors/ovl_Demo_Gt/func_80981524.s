glabel func_80981524
/* 03EB4 80981524 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 03EB8 80981528 AFA40018 */  sw      $a0, 0x0018($sp)           
/* 03EBC 8098152C AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 03EC0 80981530 AFA5001C */  sw      $a1, 0x001C($sp)           
/* 03EC4 80981534 00A02025 */  or      $a0, $a1, $zero            ## $a0 = 00000000
/* 03EC8 80981538 24050002 */  addiu   $a1, $zero, 0x0002         ## $a1 = 00000002
/* 03ECC 8098153C 0C25F9C1 */  jal     func_8097E704              
/* 03ED0 80981540 24060005 */  addiu   $a2, $zero, 0x0005         ## $a2 = 00000005
/* 03ED4 80981544 10400005 */  beq     $v0, $zero, .L8098155C     
/* 03ED8 80981548 8FA4001C */  lw      $a0, 0x001C($sp)           
/* 03EDC 8098154C 8FAF0018 */  lw      $t7, 0x0018($sp)           
/* 03EE0 80981550 240E000C */  addiu   $t6, $zero, 0x000C         ## $t6 = 0000000C
/* 03EE4 80981554 10000008 */  beq     $zero, $zero, .L80981578   
/* 03EE8 80981558 ADEE0164 */  sw      $t6, 0x0164($t7)           ## 00000164
.L8098155C:
/* 03EEC 8098155C 24050003 */  addiu   $a1, $zero, 0x0003         ## $a1 = 00000003
/* 03EF0 80981560 0C25F9C1 */  jal     func_8097E704              
/* 03EF4 80981564 24060005 */  addiu   $a2, $zero, 0x0005         ## $a2 = 00000005
/* 03EF8 80981568 10400003 */  beq     $v0, $zero, .L80981578     
/* 03EFC 8098156C 8FB90018 */  lw      $t9, 0x0018($sp)           
/* 03F00 80981570 24180011 */  addiu   $t8, $zero, 0x0011         ## $t8 = 00000011
/* 03F04 80981574 AF380164 */  sw      $t8, 0x0164($t9)           ## 00000164
.L80981578:
/* 03F08 80981578 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 03F0C 8098157C 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 03F10 80981580 03E00008 */  jr      $ra                        
/* 03F14 80981584 00000000 */  nop


