glabel func_80980504
/* 02E94 80980504 27BDFF98 */  addiu   $sp, $sp, 0xFF98           ## $sp = FFFFFF98
/* 02E98 80980508 3C0F8098 */  lui     $t7, %hi(D_80982688)       ## $t7 = 80980000
/* 02E9C 8098050C AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 02EA0 80980510 AFA40068 */  sw      $a0, 0x0068($sp)           
/* 02EA4 80980514 AFA5006C */  sw      $a1, 0x006C($sp)           
/* 02EA8 80980518 25EF2688 */  addiu   $t7, $t7, %lo(D_80982688)  ## $t7 = 80982688
/* 02EAC 8098051C 8DF90000 */  lw      $t9, 0x0000($t7)           ## 80982688
/* 02EB0 80980520 94A21D74 */  lhu     $v0, 0x1D74($a1)           ## 00001D74
/* 02EB4 80980524 27A6003C */  addiu   $a2, $sp, 0x003C           ## $a2 = FFFFFFD4
/* 02EB8 80980528 ACD90000 */  sw      $t9, 0x0000($a2)           ## FFFFFFD4
/* 02EBC 8098052C 8DF80004 */  lw      $t8, 0x0004($t7)           ## 8098268C
/* 02EC0 80980530 284102C1 */  slti    $at, $v0, 0x02C1           
/* 02EC4 80980534 3C088016 */  lui     $t0, 0x8016                ## $t0 = 80160000
/* 02EC8 80980538 ACD80004 */  sw      $t8, 0x0004($a2)           ## FFFFFFD8
/* 02ECC 8098053C 8DF90008 */  lw      $t9, 0x0008($t7)           ## 80982690
/* 02ED0 80980540 10200006 */  beq     $at, $zero, .L8098055C     
/* 02ED4 80980544 ACD90008 */  sw      $t9, 0x0008($a2)           ## FFFFFFDC
/* 02ED8 80980548 8D08FA90 */  lw      $t0, -0x0570($t0)          ## 8015FA90
/* 02EDC 8098054C 24010009 */  addiu   $at, $zero, 0x0009         ## $at = 00000009
/* 02EE0 80980550 85091456 */  lh      $t1, 0x1456($t0)           ## 80161456
/* 02EE4 80980554 5521001D */  bnel    $t1, $at, .L809805CC       
/* 02EE8 80980558 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L8098055C:
/* 02EEC 8098055C 8FA20068 */  lw      $v0, 0x0068($sp)           
/* 02EF0 80980560 3C018098 */  lui     $at, %hi(D_80982AEC)       ## $at = 80980000
/* 02EF4 80980564 C4262AEC */  lwc1    $f6, %lo(D_80982AEC)($at)  
/* 02EF8 80980568 C4440024 */  lwc1    $f4, 0x0024($v0)           ## 00000024
/* 02EFC 8098056C 3C014270 */  lui     $at, 0x4270                ## $at = 42700000
/* 02F00 80980570 44818000 */  mtc1    $at, $f16                  ## $f16 = 60.00
/* 02F04 80980574 46062200 */  add.s   $f8, $f4, $f6              
/* 02F08 80980578 3C0143C3 */  lui     $at, 0x43C3                ## $at = 43C30000
/* 02F0C 8098057C 44813000 */  mtc1    $at, $f6                   ## $f6 = 390.00
/* 02F10 80980580 240A0006 */  addiu   $t2, $zero, 0x0006         ## $t2 = 00000006
/* 02F14 80980584 E7A80050 */  swc1    $f8, 0x0050($sp)           
/* 02F18 80980588 C44A0028 */  lwc1    $f10, 0x0028($v0)          ## 00000028
/* 02F1C 8098058C 240B0002 */  addiu   $t3, $zero, 0x0002         ## $t3 = 00000002
/* 02F20 80980590 240C0023 */  addiu   $t4, $zero, 0x0023         ## $t4 = 00000023
/* 02F24 80980594 46105480 */  add.s   $f18, $f10, $f16           
/* 02F28 80980598 8FA4006C */  lw      $a0, 0x006C($sp)           
/* 02F2C 8098059C 27A50050 */  addiu   $a1, $sp, 0x0050           ## $a1 = FFFFFFE8
/* 02F30 809805A0 3C0740C0 */  lui     $a3, 0x40C0                ## $a3 = 40C00000
/* 02F34 809805A4 E7B20054 */  swc1    $f18, 0x0054($sp)          
/* 02F38 809805A8 C444002C */  lwc1    $f4, 0x002C($v0)           ## 0000002C
/* 02F3C 809805AC AFAC0018 */  sw      $t4, 0x0018($sp)           
/* 02F40 809805B0 AFAB0014 */  sw      $t3, 0x0014($sp)           
/* 02F44 809805B4 46062200 */  add.s   $f8, $f4, $f6              
/* 02F48 809805B8 AFAA0010 */  sw      $t2, 0x0010($sp)           
/* 02F4C 809805BC 24420024 */  addiu   $v0, $v0, 0x0024           ## $v0 = 00000024
/* 02F50 809805C0 0C25F5F6 */  jal     func_8097D7D8              
/* 02F54 809805C4 E7A80058 */  swc1    $f8, 0x0058($sp)           
/* 02F58 809805C8 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L809805CC:
/* 02F5C 809805CC 27BD0068 */  addiu   $sp, $sp, 0x0068           ## $sp = 00000000
/* 02F60 809805D0 03E00008 */  jr      $ra                        
/* 02F64 809805D4 00000000 */  nop


