.late_rodata
glabel D_80937A10
 .word 0x3F555555

.text
glabel func_80933B84
/* 075B4 80933B84 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 075B8 80933B88 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 075BC 80933B8C 0C03F66B */  jal     Math_Rand_ZeroOne
              ## Rand.Next() float
/* 075C0 80933B90 AFA40018 */  sw      $a0, 0x0018($sp)           
/* 075C4 80933B94 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 075C8 80933B98 3C0140C0 */  lui     $at, 0x40C0                ## $at = 40C00000
/* 075CC 80933B9C 44812000 */  mtc1    $at, $f4                   ## $f4 = 6.00
/* 075D0 80933BA0 8C8E011C */  lw      $t6, 0x011C($a0)           ## 0000011C
/* 075D4 80933BA4 3C198093 */  lui     $t9, %hi(D_8093746C)       ## $t9 = 80930000
/* 075D8 80933BA8 46040302 */  mul.s   $f12, $f0, $f4             
/* 075DC 80933BAC 85CF001C */  lh      $t7, 0x001C($t6)           ## 0000001C
/* 075E0 80933BB0 24010008 */  addiu   $at, $zero, 0x0008         ## $at = 00000008
/* 075E4 80933BB4 000FC080 */  sll     $t8, $t7,  2               
/* 075E8 80933BB8 0338C821 */  addu    $t9, $t9, $t8              
/* 075EC 80933BBC 8F39746C */  lw      $t9, %lo(D_8093746C)($t9)  
/* 075F0 80933BC0 46006086 */  mov.s   $f2, $f12                  
/* 075F4 80933BC4 1721000A */  bne     $t9, $at, .L80933BF0       
/* 075F8 80933BC8 3C014080 */  lui     $at, 0x4080                ## $at = 40800000
/* 075FC 80933BCC 44810000 */  mtc1    $at, $f0                   ## $f0 = 4.00
/* 07600 80933BD0 3C018093 */  lui     $at, %hi(D_80937A10)       ## $at = 80930000
/* 07604 80933BD4 C4267A10 */  lwc1    $f6, %lo(D_80937A10)($at)  
/* 07608 80933BD8 46066082 */  mul.s   $f2, $f12, $f6             
/* 0760C 80933BDC 4602003C */  c.lt.s  $f0, $f2                   
/* 07610 80933BE0 00000000 */  nop
/* 07614 80933BE4 45020003 */  bc1fl   .L80933BF4                 
/* 07618 80933BE8 4600120D */  trunc.w.s $f8, $f2                   
/* 0761C 80933BEC 46000086 */  mov.s   $f2, $f0                   
.L80933BF0:
/* 07620 80933BF0 4600120D */  trunc.w.s $f8, $f2                   
.L80933BF4:
/* 07624 80933BF4 24010001 */  addiu   $at, $zero, 0x0001         ## $at = 00000001
/* 07628 80933BF8 44024000 */  mfc1    $v0, $f8                   
/* 0762C 80933BFC 00000000 */  nop
/* 07630 80933C00 14400005 */  bne     $v0, $zero, .L80933C18     
/* 07634 80933C04 00000000 */  nop
/* 07638 80933C08 0C24C1D2 */  jal     func_80930748              
/* 0763C 80933C0C 00000000 */  nop
/* 07640 80933C10 10000018 */  beq     $zero, $zero, .L80933C74   
/* 07644 80933C14 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80933C18:
/* 07648 80933C18 54410006 */  bnel    $v0, $at, .L80933C34       
/* 0764C 80933C1C 24010002 */  addiu   $at, $zero, 0x0002         ## $at = 00000002
/* 07650 80933C20 0C24C2C6 */  jal     func_80930B18              
/* 07654 80933C24 00000000 */  nop
/* 07658 80933C28 10000012 */  beq     $zero, $zero, .L80933C74   
/* 0765C 80933C2C 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 07660 80933C30 24010002 */  addiu   $at, $zero, 0x0002         ## $at = 00000002
.L80933C34:
/* 07664 80933C34 54410006 */  bnel    $v0, $at, .L80933C50       
/* 07668 80933C38 24010005 */  addiu   $at, $zero, 0x0005         ## $at = 00000005
/* 0766C 80933C3C 0C24C3E0 */  jal     func_80930F80              
/* 07670 80933C40 00000000 */  nop
/* 07674 80933C44 1000000B */  beq     $zero, $zero, .L80933C74   
/* 07678 80933C48 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 0767C 80933C4C 24010005 */  addiu   $at, $zero, 0x0005         ## $at = 00000005
.L80933C50:
/* 07680 80933C50 14410005 */  bne     $v0, $at, .L80933C68       
/* 07684 80933C54 00000000 */  nop
/* 07688 80933C58 0C24C484 */  jal     func_80931210              
/* 0768C 80933C5C 00000000 */  nop
/* 07690 80933C60 10000004 */  beq     $zero, $zero, .L80933C74   
/* 07694 80933C64 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80933C68:
/* 07698 80933C68 0C24C628 */  jal     func_809318A0              
/* 0769C 80933C6C 00000000 */  nop
/* 076A0 80933C70 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80933C74:
/* 076A4 80933C74 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 076A8 80933C78 03E00008 */  jr      $ra                        
/* 076AC 80933C7C 00000000 */  nop
