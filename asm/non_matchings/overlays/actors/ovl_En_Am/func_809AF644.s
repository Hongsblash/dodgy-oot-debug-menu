glabel func_809AF644
/* 01724 809AF644 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 01728 809AF648 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 0172C 809AF64C 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 01730 809AF650 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 01734 809AF654 AFA5002C */  sw      $a1, 0x002C($sp)           
/* 01738 809AF658 86050032 */  lh      $a1, 0x0032($s0)           ## 00000032
/* 0173C 809AF65C AFA00010 */  sw      $zero, 0x0010($sp)         
/* 01740 809AF660 248400B6 */  addiu   $a0, $a0, 0x00B6           ## $a0 = 000000B6
/* 01744 809AF664 24060001 */  addiu   $a2, $zero, 0x0001         ## $a2 = 00000001
/* 01748 809AF668 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS
              
/* 0174C 809AF66C 24070FA0 */  addiu   $a3, $zero, 0x0FA0         ## $a3 = 00000FA0
/* 01750 809AF670 44801000 */  mtc1    $zero, $f2                 ## $f2 = 0.00
/* 01754 809AF674 C6000068 */  lwc1    $f0, 0x0068($s0)           ## 00000068
/* 01758 809AF678 3C013F00 */  lui     $at, 0x3F00                ## $at = 3F000000
/* 0175C 809AF67C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 01760 809AF680 4602003C */  c.lt.s  $f0, $f2                   
/* 01764 809AF684 3C06C110 */  lui     $a2, 0xC110                ## $a2 = C1100000
/* 01768 809AF688 45020006 */  bc1fl   .L809AF6A4                 
/* 0176C 809AF68C C6080060 */  lwc1    $f8, 0x0060($s0)           ## 00000060
/* 01770 809AF690 44812000 */  mtc1    $at, $f4                   ## $f4 = 0.50
/* 01774 809AF694 00000000 */  nop
/* 01778 809AF698 46040180 */  add.s   $f6, $f0, $f4              
/* 0177C 809AF69C E6060068 */  swc1    $f6, 0x0068($s0)           ## 00000068
/* 01780 809AF6A0 C6080060 */  lwc1    $f8, 0x0060($s0)           ## 00000060
.L809AF6A4:
/* 01784 809AF6A4 8FA5002C */  lw      $a1, 0x002C($sp)           
/* 01788 809AF6A8 4602403E */  c.le.s  $f8, $f2                   
/* 0178C 809AF6AC 00000000 */  nop
/* 01790 809AF6B0 45020008 */  bc1fl   .L809AF6D4                 
/* 01794 809AF6B4 920E0114 */  lbu     $t6, 0x0114($s0)           ## 00000114
/* 01798 809AF6B8 0C26B7CA */  jal     func_809ADF28              
/* 0179C 809AF6BC 86070032 */  lh      $a3, 0x0032($s0)           ## 00000032
/* 017A0 809AF6C0 44801000 */  mtc1    $zero, $f2                 ## $f2 = 0.00
/* 017A4 809AF6C4 54400003 */  bnel    $v0, $zero, .L809AF6D4     
/* 017A8 809AF6C8 920E0114 */  lbu     $t6, 0x0114($s0)           ## 00000114
/* 017AC 809AF6CC E6020068 */  swc1    $f2, 0x0068($s0)           ## 00000068
/* 017B0 809AF6D0 920E0114 */  lbu     $t6, 0x0114($s0)           ## 00000114
.L809AF6D4:
/* 017B4 809AF6D4 55C0000C */  bnel    $t6, $zero, .L809AF708     
/* 017B8 809AF6D8 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 017BC 809AF6DC 920F00AF */  lbu     $t7, 0x00AF($s0)           ## 000000AF
/* 017C0 809AF6E0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 017C4 809AF6E4 11E00005 */  beq     $t7, $zero, .L809AF6FC     
/* 017C8 809AF6E8 00000000 */  nop
/* 017CC 809AF6EC 0C26B94E */  jal     func_809AE538              
/* 017D0 809AF6F0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 017D4 809AF6F4 10000004 */  beq     $zero, $zero, .L809AF708   
/* 017D8 809AF6F8 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L809AF6FC:
/* 017DC 809AF6FC 0C26B9C7 */  jal     func_809AE71C              
/* 017E0 809AF700 8FA5002C */  lw      $a1, 0x002C($sp)           
/* 017E4 809AF704 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L809AF708:
/* 017E8 809AF708 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 017EC 809AF70C 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 017F0 809AF710 03E00008 */  jr      $ra                        
/* 017F4 809AF714 00000000 */  nop


