glabel func_80A891C4
/* 003B4 80A891C4 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 003B8 80A891C8 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 003BC 80A891CC AFA40018 */  sw      $a0, 0x0018($sp)           
/* 003C0 80A891D0 24A420D8 */  addiu   $a0, $a1, 0x20D8           ## $a0 = 000020D8
/* 003C4 80A891D4 0C042F6F */  jal     func_8010BDBC              
/* 003C8 80A891D8 AFA5001C */  sw      $a1, 0x001C($sp)           
/* 003CC 80A891DC 24010004 */  addiu   $at, $zero, 0x0004         ## $at = 00000004
/* 003D0 80A891E0 14410028 */  bne     $v0, $at, .L80A89284       
/* 003D4 80A891E4 8FA6001C */  lw      $a2, 0x001C($sp)           
/* 003D8 80A891E8 00C02025 */  or      $a0, $a2, $zero            ## $a0 = 00000000
/* 003DC 80A891EC 0C041AF2 */  jal     func_80106BC8              
/* 003E0 80A891F0 AFA6001C */  sw      $a2, 0x001C($sp)           
/* 003E4 80A891F4 10400023 */  beq     $v0, $zero, .L80A89284     
/* 003E8 80A891F8 8FA6001C */  lw      $a2, 0x001C($sp)           
/* 003EC 80A891FC 3C020001 */  lui     $v0, 0x0001                ## $v0 = 00010000
/* 003F0 80A89200 00461021 */  addu    $v0, $v0, $a2              
/* 003F4 80A89204 904204BD */  lbu     $v0, 0x04BD($v0)           ## 000104BD
/* 003F8 80A89208 24010001 */  addiu   $at, $zero, 0x0001         ## $at = 00000001
/* 003FC 80A8920C 3C0E8016 */  lui     $t6, 0x8016                ## $t6 = 80160000
/* 00400 80A89210 10400005 */  beq     $v0, $zero, .L80A89228     
/* 00404 80A89214 00000000 */  nop
/* 00408 80A89218 10410016 */  beq     $v0, $at, .L80A89274       
/* 0040C 80A8921C 00C02025 */  or      $a0, $a2, $zero            ## $a0 = 00000000
/* 00410 80A89220 10000019 */  beq     $zero, $zero, .L80A89288   
/* 00414 80A89224 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80A89228:
/* 00418 80A89228 85CEE694 */  lh      $t6, -0x196C($t6)          ## 8015E694
/* 0041C 80A8922C 00C02025 */  or      $a0, $a2, $zero            ## $a0 = 00000000
/* 00420 80A89230 29C100C8 */  slti    $at, $t6, 0x00C8           
/* 00424 80A89234 10200007 */  beq     $at, $zero, .L80A89254     
/* 00428 80A89238 00000000 */  nop
/* 0042C 80A8923C 0C042DC8 */  jal     func_8010B720              
/* 00430 80A89240 24056075 */  addiu   $a1, $zero, 0x6075         ## $a1 = 00006075
/* 00434 80A89244 0C2A2402 */  jal     func_80A89008              
/* 00438 80A89248 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 0043C 80A8924C 1000000E */  beq     $zero, $zero, .L80A89288   
/* 00440 80A89250 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80A89254:
/* 00444 80A89254 0C021CC3 */  jal     Rupees_ChangeBy              
/* 00448 80A89258 2404FF38 */  addiu   $a0, $zero, 0xFF38         ## $a0 = FFFFFF38
/* 0044C 80A8925C 3C0580A9 */  lui     $a1, %hi(func_80A89160)    ## $a1 = 80A90000
/* 00450 80A89260 24A59160 */  addiu   $a1, $a1, %lo(func_80A89160) ## $a1 = 80A89160
/* 00454 80A89264 0C2A2384 */  jal     func_80A88E10              
/* 00458 80A89268 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 0045C 80A8926C 10000006 */  beq     $zero, $zero, .L80A89288   
/* 00460 80A89270 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80A89274:
/* 00464 80A89274 0C042DC8 */  jal     func_8010B720              
/* 00468 80A89278 24056074 */  addiu   $a1, $zero, 0x6074         ## $a1 = 00006074
/* 0046C 80A8927C 0C2A2402 */  jal     func_80A89008              
/* 00470 80A89280 8FA40018 */  lw      $a0, 0x0018($sp)           
.L80A89284:
/* 00474 80A89284 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80A89288:
/* 00478 80A89288 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 0047C 80A8928C 03E00008 */  jr      $ra                        
/* 00480 80A89290 00000000 */  nop


