glabel func_80A5760C
/* 0040C 80A5760C 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 00410 80A57610 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 00414 80A57614 00803825 */  or      $a3, $a0, $zero            ## $a3 = 00000000
/* 00418 80A57618 3C050600 */  lui     $a1, 0x0600                ## $a1 = 06000000
/* 0041C 80A5761C 24A526C4 */  addiu   $a1, $a1, 0x26C4           ## $a1 = 060026C4
/* 00420 80A57620 AFA70018 */  sw      $a3, 0x0018($sp)           
/* 00424 80A57624 2484014C */  addiu   $a0, $a0, 0x014C           ## $a0 = 0000014C
/* 00428 80A57628 0C029490 */  jal     func_800A5240              
/* 0042C 80A5762C 3C06C040 */  lui     $a2, 0xC040                ## $a2 = C0400000
/* 00430 80A57630 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 00434 80A57634 240E0025 */  addiu   $t6, $zero, 0x0025         ## $t6 = 00000025
/* 00438 80A57638 24053880 */  addiu   $a1, $zero, 0x3880         ## $a1 = 00003880
/* 0043C 80A5763C 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 00440 80A57640 A48E0256 */  sh      $t6, 0x0256($a0)           ## 00000256
/* 00444 80A57644 8FA70018 */  lw      $a3, 0x0018($sp)           
/* 00448 80A57648 3C0D80A5 */  lui     $t5, %hi(func_80A57EF8)    ## $t5 = 80A50000
/* 0044C 80A5764C 25AD7EF8 */  addiu   $t5, $t5, %lo(func_80A57EF8) ## $t5 = 80A57EF8
/* 00450 80A57650 90EF0225 */  lbu     $t7, 0x0225($a3)           ## 00000225
/* 00454 80A57654 84E3001C */  lh      $v1, 0x001C($a3)           ## 0000001C
/* 00458 80A57658 31F8FFFE */  andi    $t8, $t7, 0xFFFE           ## $t8 = 00000000
/* 0045C 80A5765C 18600022 */  blez    $v1, .L80A576E8            
/* 00460 80A57660 A0F80225 */  sb      $t8, 0x0225($a3)           ## 00000225
/* 00464 80A57664 28610004 */  slti    $at, $v1, 0x0004           
/* 00468 80A57668 50200020 */  beql    $at, $zero, .L80A576EC     
/* 0046C 80A5766C ACED0190 */  sw      $t5, 0x0190($a3)           ## 00000190
/* 00470 80A57670 90F90002 */  lbu     $t9, 0x0002($a3)           ## 00000002
/* 00474 80A57674 24010005 */  addiu   $at, $zero, 0x0005         ## $at = 00000005
/* 00478 80A57678 3C0580A6 */  lui     $a1, %hi(D_80A58A34)       ## $a1 = 80A60000
/* 0047C 80A5767C 1721001A */  bne     $t9, $at, .L80A576E8       
/* 00480 80A57680 24A58A34 */  addiu   $a1, $a1, %lo(D_80A58A34)  ## $a1 = 80A58A34
/* 00484 80A57684 84A20000 */  lh      $v0, 0x0000($a1)           ## 80A58A34
/* 00488 80A57688 2401FFFC */  addiu   $at, $zero, 0xFFFC         ## $at = FFFFFFFC
/* 0048C 80A5768C 3C0C80A5 */  lui     $t4, %hi(func_80A57F5C)    ## $t4 = 80A50000
/* 00490 80A57690 54410005 */  bnel    $v0, $at, .L80A576A8       
/* 00494 80A57694 24440001 */  addiu   $a0, $v0, 0x0001           ## $a0 = 00000001
/* 00498 80A57698 A4A00000 */  sh      $zero, 0x0000($a1)         ## 80A58A34
/* 0049C 80A5769C 84A20000 */  lh      $v0, 0x0000($a1)           ## 80A58A34
/* 004A0 80A576A0 84E3001C */  lh      $v1, 0x001C($a3)           ## 0000001C
/* 004A4 80A576A4 24440001 */  addiu   $a0, $v0, 0x0001           ## $a0 = 00000001
.L80A576A8:
/* 004A8 80A576A8 14830003 */  bne     $a0, $v1, .L80A576B8       
/* 004AC 80A576AC 00000000 */  nop
/* 004B0 80A576B0 10000007 */  beq     $zero, $zero, .L80A576D0   
/* 004B4 80A576B4 A4A40000 */  sh      $a0, 0x0000($a1)           ## 80A58A34
.L80A576B8:
/* 004B8 80A576B8 18400003 */  blez    $v0, .L80A576C8            
/* 004BC 80A576BC 00024023 */  subu    $t0, $zero, $v0            
/* 004C0 80A576C0 A4A80000 */  sh      $t0, 0x0000($a1)           ## 80A58A34
/* 004C4 80A576C4 84A20000 */  lh      $v0, 0x0000($a1)           ## 80A58A34
.L80A576C8:
/* 004C8 80A576C8 2449FFFF */  addiu   $t1, $v0, 0xFFFF           ## $t1 = FFFFFFFF
/* 004CC 80A576CC A4A90000 */  sh      $t1, 0x0000($a1)           ## 80A58A34
.L80A576D0:
/* 004D0 80A576D0 8CEA0004 */  lw      $t2, 0x0004($a3)           ## 00000004
/* 004D4 80A576D4 258C7F5C */  addiu   $t4, $t4, %lo(func_80A57F5C) ## $t4 = 80A57F5C
/* 004D8 80A576D8 ACEC0190 */  sw      $t4, 0x0190($a3)           ## 00000190
/* 004DC 80A576DC 354B0010 */  ori     $t3, $t2, 0x0010           ## $t3 = 00000010
/* 004E0 80A576E0 10000002 */  beq     $zero, $zero, .L80A576EC   
/* 004E4 80A576E4 ACEB0004 */  sw      $t3, 0x0004($a3)           ## 00000004
.L80A576E8:
/* 004E8 80A576E8 ACED0190 */  sw      $t5, 0x0190($a3)           ## 00000190
.L80A576EC:
/* 004EC 80A576EC 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 004F0 80A576F0 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 004F4 80A576F4 03E00008 */  jr      $ra                        
/* 004F8 80A576F8 00000000 */  nop


