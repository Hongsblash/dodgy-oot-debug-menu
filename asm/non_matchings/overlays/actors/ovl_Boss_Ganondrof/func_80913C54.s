glabel func_80913C54
/* 03614 80913C54 27BDFFD0 */  addiu   $sp, $sp, 0xFFD0           ## $sp = FFFFFFD0
/* 03618 80913C58 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 0361C 80913C5C AFB00018 */  sw      $s0, 0x0018($sp)           
/* 03620 80913C60 AFA50034 */  sw      $a1, 0x0034($sp)           
/* 03624 80913C64 848201A0 */  lh      $v0, 0x01A0($a0)           ## 000001A0
/* 03628 80913C68 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 0362C 80913C6C 8C87011C */  lw      $a3, 0x011C($a0)           ## 0000011C
/* 03630 80913C70 10400007 */  beq     $v0, $zero, .L80913C90     
/* 03634 80913C74 244EFFFF */  addiu   $t6, $v0, 0xFFFF           ## $t6 = FFFFFFFF
/* 03638 80913C78 908F04F1 */  lbu     $t7, 0x04F1($a0)           ## 000004F1
/* 0363C 80913C7C A48E01A0 */  sh      $t6, 0x01A0($a0)           ## 000001A0
/* 03640 80913C80 A08001C7 */  sb      $zero, 0x01C7($a0)         ## 000001C7
/* 03644 80913C84 31F8FFFD */  andi    $t8, $t7, 0xFFFD           ## $t8 = 00000000
/* 03648 80913C88 1000007B */  beq     $zero, $zero, .L80913E78   
/* 0364C 80913C8C A09804F1 */  sb      $t8, 0x04F1($a0)           ## 000004F1
.L80913C90:
/* 03650 80913C90 920304F1 */  lbu     $v1, 0x04F1($s0)           ## 000004F1
/* 03654 80913C94 30620002 */  andi    $v0, $v1, 0x0002           ## $v0 = 00000000
/* 03658 80913C98 50400005 */  beql    $v0, $zero, .L80913CB0     
/* 0365C 80913C9C 920801C7 */  lbu     $t0, 0x01C7($s0)           ## 000001C7
/* 03660 80913CA0 821900AF */  lb      $t9, 0x00AF($s0)           ## 000000AF
/* 03664 80913CA4 1F200004 */  bgtz    $t9, .L80913CB8            
/* 03668 80913CA8 00000000 */  nop
/* 0366C 80913CAC 920801C7 */  lbu     $t0, 0x01C7($s0)           ## 000001C7
.L80913CB0:
/* 03670 80913CB0 51000072 */  beql    $t0, $zero, .L80913E7C     
/* 03674 80913CB4 8FBF001C */  lw      $ra, 0x001C($sp)           
.L80913CB8:
/* 03678 80913CB8 50400006 */  beql    $v0, $zero, .L80913CD4     
/* 0367C 80913CBC 920A01C9 */  lbu     $t2, 0x01C9($s0)           ## 000001C9
/* 03680 80913CC0 8E06051C */  lw      $a2, 0x051C($s0)           ## 0000051C
/* 03684 80913CC4 3069FFFD */  andi    $t1, $v1, 0xFFFD           ## $t1 = 00000000
/* 03688 80913CC8 A20904F1 */  sb      $t1, 0x04F1($s0)           ## 000004F1
/* 0368C 80913CCC AFA60024 */  sw      $a2, 0x0024($sp)           
/* 03690 80913CD0 920A01C9 */  lbu     $t2, 0x01C9($s0)           ## 000001C9
.L80913CD4:
/* 03694 80913CD4 8FA60024 */  lw      $a2, 0x0024($sp)           
/* 03698 80913CD8 11400055 */  beq     $t2, $zero, .L80913E30     
/* 0369C 80913CDC 00000000 */  nop
/* 036A0 80913CE0 50400014 */  beql    $v0, $zero, .L80913D34     
/* 036A4 80913CE4 8E180190 */  lw      $t8, 0x0190($s0)           ## 00000190
/* 036A8 80913CE8 8E0C0190 */  lw      $t4, 0x0190($s0)           ## 00000190
/* 036AC 80913CEC 3C0B8091 */  lui     $t3, %hi(func_809122A4)    ## $t3 = 80910000
/* 036B0 80913CF0 256B22A4 */  addiu   $t3, $t3, %lo(func_809122A4) ## $t3 = 809122A4
/* 036B4 80913CF4 516C000F */  beql    $t3, $t4, .L80913D34       
/* 036B8 80913CF8 8E180190 */  lw      $t8, 0x0190($s0)           ## 00000190
/* 036BC 80913CFC 8CCD0000 */  lw      $t5, 0x0000($a2)           ## 00000000
/* 036C0 80913D00 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 036C4 80913D04 3421F8A4 */  ori     $at, $at, 0xF8A4           ## $at = 0001F8A4
/* 036C8 80913D08 01A17024 */  and     $t6, $t5, $at              
/* 036CC 80913D0C 11C00008 */  beq     $t6, $zero, .L80913D30     
/* 036D0 80913D10 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 036D4 80913D14 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 036D8 80913D18 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 036DC 80913D1C 3C048091 */  lui     $a0, %hi(D_80915308)       ## $a0 = 80910000
/* 036E0 80913D20 0C00084C */  jal     osSyncPrintf
              
/* 036E4 80913D24 24845308 */  addiu   $a0, $a0, %lo(D_80915308)  ## $a0 = 80915308
/* 036E8 80913D28 10000053 */  beq     $zero, $zero, .L80913E78   
/* 036EC 80913D2C A20001C7 */  sb      $zero, 0x01C7($s0)         ## 000001C7
.L80913D30:
/* 036F0 80913D30 8E180190 */  lw      $t8, 0x0190($s0)           ## 00000190
.L80913D34:
/* 036F4 80913D34 3C0F8091 */  lui     $t7, %hi(func_80912594)    ## $t7 = 80910000
/* 036F8 80913D38 25EF2594 */  addiu   $t7, $t7, %lo(func_80912594) ## $t7 = 80912594
/* 036FC 80913D3C 11F80038 */  beq     $t7, $t8, .L80913E20       
/* 03700 80913D40 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 03704 80913D44 921901C7 */  lbu     $t9, 0x01C7($s0)           ## 000001C7
/* 03708 80913D48 57200023 */  bnel    $t9, $zero, .L80913DD8     
/* 0370C 80913D4C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 03710 80913D50 8CC40000 */  lw      $a0, 0x0000($a2)           ## 00000000
/* 03714 80913D54 30880080 */  andi    $t0, $a0, 0x0080           ## $t0 = 00000000
/* 03718 80913D58 55000048 */  bnel    $t0, $zero, .L80913E7C     
/* 0371C 80913D5C 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 03720 80913D60 A3A00022 */  sb      $zero, 0x0022($sp)         
/* 03724 80913D64 0C018D74 */  jal     func_800635D0              
/* 03728 80913D68 AFA70028 */  sw      $a3, 0x0028($sp)           
/* 0372C 80913D6C 93A50022 */  lbu     $a1, 0x0022($sp)           
/* 03730 80913D70 8FA70028 */  lw      $a3, 0x0028($sp)           
/* 03734 80913D74 14400003 */  bne     $v0, $zero, .L80913D84     
/* 03738 80913D78 304400FF */  andi    $a0, $v0, 0x00FF           ## $a0 = 00000000
/* 0373C 80913D7C 10000002 */  beq     $zero, $zero, .L80913D88   
/* 03740 80913D80 24040002 */  addiu   $a0, $zero, 0x0002         ## $a0 = 00000002
.L80913D84:
/* 03744 80913D84 24050001 */  addiu   $a1, $zero, 0x0001         ## $a1 = 00000001
.L80913D88:
/* 03748 80913D88 920300AF */  lbu     $v1, 0x00AF($s0)           ## 000000AF
/* 0374C 80913D8C 00031600 */  sll     $v0, $v1, 24               
/* 03750 80913D90 00021603 */  sra     $v0, $v0, 24               
/* 03754 80913D94 28410003 */  slti    $at, $v0, 0x0003           
/* 03758 80913D98 10200003 */  beq     $at, $zero, .L80913DA8     
/* 0375C 80913D9C 00644823 */  subu    $t1, $v1, $a0              
/* 03760 80913DA0 10A00003 */  beq     $a1, $zero, .L80913DB0     
/* 03764 80913DA4 00000000 */  nop
.L80913DA8:
/* 03768 80913DA8 A20900AF */  sb      $t1, 0x00AF($s0)           ## 000000AF
/* 0376C 80913DAC 820200AF */  lb      $v0, 0x00AF($s0)           ## 000000AF
.L80913DB0:
/* 03770 80913DB0 1C400008 */  bgtz    $v0, .L80913DD4            
/* 03774 80913DB4 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 03778 80913DB8 0C244B25 */  jal     func_80912C94              
/* 0377C 80913DBC 8FA50034 */  lw      $a1, 0x0034($sp)           
/* 03780 80913DC0 8FA40034 */  lw      $a0, 0x0034($sp)           
/* 03784 80913DC4 0C00CB1F */  jal     func_80032C7C              
/* 03788 80913DC8 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 0378C 80913DCC 1000002B */  beq     $zero, $zero, .L80913E7C   
/* 03790 80913DD0 8FBF001C */  lw      $ra, 0x001C($sp)           
.L80913DD4:
/* 03794 80913DD4 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
.L80913DD8:
/* 03798 80913DD8 8FA50034 */  lw      $a1, 0x0034($sp)           
/* 0379C 80913DDC 0C244878 */  jal     func_809121E0              
/* 037A0 80913DE0 AFA70028 */  sw      $a3, 0x0028($sp)           
/* 037A4 80913DE4 920A01C7 */  lbu     $t2, 0x01C7($s0)           ## 000001C7
/* 037A8 80913DE8 8FA70028 */  lw      $a3, 0x0028($sp)           
/* 037AC 80913DEC 240C000A */  addiu   $t4, $zero, 0x000A         ## $t4 = 0000000A
/* 037B0 80913DF0 29410002 */  slti    $at, $t2, 0x0002           
/* 037B4 80913DF4 14200003 */  bne     $at, $zero, .L80913E04     
/* 037B8 80913DF8 240D0014 */  addiu   $t5, $zero, 0x0014         ## $t5 = 00000014
/* 037BC 80913DFC 240B0078 */  addiu   $t3, $zero, 0x0078         ## $t3 = 00000078
/* 037C0 80913E00 A60B01BC */  sh      $t3, 0x01BC($s0)           ## 000001BC
.L80913E04:
/* 037C4 80913E04 A60C01A0 */  sh      $t4, 0x01A0($s0)           ## 000001A0
/* 037C8 80913E08 A4ED01DE */  sh      $t5, 0x01DE($a3)           ## 000001DE
/* 037CC 80913E0C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 037D0 80913E10 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 037D4 80913E14 240538AE */  addiu   $a1, $zero, 0x38AE         ## $a1 = 000038AE
/* 037D8 80913E18 10000017 */  beq     $zero, $zero, .L80913E78   
/* 037DC 80913E1C A20001C7 */  sb      $zero, 0x01C7($s0)         ## 000001C7
.L80913E20:
/* 037E0 80913E20 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 037E4 80913E24 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 037E8 80913E28 10000013 */  beq     $zero, $zero, .L80913E78   
/* 037EC 80913E2C A20001C7 */  sb      $zero, 0x01C7($s0)         ## 000001C7
.L80913E30:
/* 037F0 80913E30 50400011 */  beql    $v0, $zero, .L80913E78     
/* 037F4 80913E34 A20001C7 */  sb      $zero, 0x01C7($s0)         ## 000001C7
/* 037F8 80913E38 8CCE0000 */  lw      $t6, 0x0000($a2)           ## 00000000
/* 037FC 80913E3C 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 03800 80913E40 3421F8A4 */  ori     $at, $at, 0xF8A4           ## $at = 0001F8A4
/* 03804 80913E44 01C17824 */  and     $t7, $t6, $at              
/* 03808 80913E48 11E0000A */  beq     $t7, $zero, .L80913E74     
/* 0380C 80913E4C 2418000A */  addiu   $t8, $zero, 0x000A         ## $t8 = 0000000A
/* 03810 80913E50 921900AF */  lbu     $t9, 0x00AF($s0)           ## 000000AF
/* 03814 80913E54 A61801A0 */  sh      $t8, 0x01A0($s0)           ## 000001A0
/* 03818 80913E58 24090014 */  addiu   $t1, $zero, 0x0014         ## $t1 = 00000014
/* 0381C 80913E5C 2728FFFE */  addiu   $t0, $t9, 0xFFFE           ## $t0 = FFFFFFFE
/* 03820 80913E60 A20800AF */  sb      $t0, 0x00AF($s0)           ## 000000AF
/* 03824 80913E64 A4E901DE */  sh      $t1, 0x01DE($a3)           ## 000001DE
/* 03828 80913E68 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 0382C 80913E6C 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 03830 80913E70 240538AE */  addiu   $a1, $zero, 0x38AE         ## $a1 = 000038AE
.L80913E74:
/* 03834 80913E74 A20001C7 */  sb      $zero, 0x01C7($s0)         ## 000001C7
.L80913E78:
/* 03838 80913E78 8FBF001C */  lw      $ra, 0x001C($sp)           
.L80913E7C:
/* 0383C 80913E7C 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 03840 80913E80 27BD0030 */  addiu   $sp, $sp, 0x0030           ## $sp = 00000000
/* 03844 80913E84 03E00008 */  jr      $ra                        
/* 03848 80913E88 00000000 */  nop


