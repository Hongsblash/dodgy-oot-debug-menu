glabel EnKo_Draw
/* 02F50 80A99D00 27BDFF98 */  addiu   $sp, $sp, 0xFF98           ## $sp = FFFFFF98
/* 02F54 80A99D04 AFBF002C */  sw      $ra, 0x002C($sp)           
/* 02F58 80A99D08 AFB20028 */  sw      $s2, 0x0028($sp)           
/* 02F5C 80A99D0C AFB10024 */  sw      $s1, 0x0024($sp)           
/* 02F60 80A99D10 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 02F64 80A99D14 848F001C */  lh      $t7, 0x001C($a0)           ## 0000001C
/* 02F68 80A99D18 2403000B */  addiu   $v1, $zero, 0x000B         ## $v1 = 0000000B
/* 02F6C 80A99D1C 3C0280AA */  lui     $v0, %hi(D_80A9A500)       ## $v0 = 80AA0000
/* 02F70 80A99D20 31F800FF */  andi    $t8, $t7, 0x00FF           ## $t8 = 00000000
/* 02F74 80A99D24 03030019 */  multu   $t8, $v1                   
/* 02F78 80A99D28 2442A500 */  addiu   $v0, $v0, %lo(D_80A9A500)  ## $v0 = 80A9A500
/* 02F7C 80A99D2C 27AE0060 */  addiu   $t6, $sp, 0x0060           ## $t6 = FFFFFFF8
/* 02F80 80A99D30 27AB005C */  addiu   $t3, $sp, 0x005C           ## $t3 = FFFFFFF4
/* 02F84 80A99D34 24090001 */  addiu   $t1, $zero, 0x0001         ## $t1 = 00000001
/* 02F88 80A99D38 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 02F8C 80A99D3C 3C0680AA */  lui     $a2, %hi(D_80A9A79C)       ## $a2 = 80AA0000
/* 02F90 80A99D40 00A09025 */  or      $s2, $a1, $zero            ## $s2 = 00000000
/* 02F94 80A99D44 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 02F98 80A99D48 24C6A79C */  addiu   $a2, $a2, %lo(D_80A9A79C)  ## $a2 = 80A9A79C
/* 02F9C 80A99D4C 0000C812 */  mflo    $t9                        
/* 02FA0 80A99D50 00594021 */  addu    $t0, $v0, $t9              
/* 02FA4 80A99D54 890A0002 */  lwl     $t2, 0x0002($t0)           ## 00000002
/* 02FA8 80A99D58 990A0005 */  lwr     $t2, 0x0005($t0)           ## 00000005
/* 02FAC 80A99D5C 4448F800 */  cfc1    $t0, $31
/* 02FB0 80A99D60 44C9F800 */  ctc1    $t1, $31
/* 02FB4 80A99D64 ADCA0000 */  sw      $t2, 0x0000($t6)           ## FFFFFFF8
/* 02FB8 80A99D68 848C001C */  lh      $t4, 0x001C($a0)           ## 0000001C
/* 02FBC 80A99D6C 2407082F */  addiu   $a3, $zero, 0x082F         ## $a3 = 0000082F
/* 02FC0 80A99D70 318D00FF */  andi    $t5, $t4, 0x00FF           ## $t5 = 00000000
/* 02FC4 80A99D74 01A30019 */  multu   $t5, $v1                   
/* 02FC8 80A99D78 00007812 */  mflo    $t7                        
/* 02FCC 80A99D7C 004FC021 */  addu    $t8, $v0, $t7              
/* 02FD0 80A99D80 8B0E0007 */  lwl     $t6, 0x0007($t8)           ## 00000007
/* 02FD4 80A99D84 9B0E000A */  lwr     $t6, 0x000A($t8)           ## 0000000A
/* 02FD8 80A99D88 AD6E0000 */  sw      $t6, 0x0000($t3)           ## FFFFFFF4
/* 02FDC 80A99D8C C4840220 */  lwc1    $f4, 0x0220($a0)           ## 00000220
/* 02FE0 80A99D90 27A40048 */  addiu   $a0, $sp, 0x0048           ## $a0 = FFFFFFE0
/* 02FE4 80A99D94 460021A4 */  cvt.w.s $f6, $f4                   
/* 02FE8 80A99D98 4449F800 */  cfc1    $t1, $31
/* 02FEC 80A99D9C 00000000 */  nop
/* 02FF0 80A99DA0 31290078 */  andi    $t1, $t1, 0x0078           ## $t1 = 00000000
/* 02FF4 80A99DA4 51200013 */  beql    $t1, $zero, .L80A99DF4     
/* 02FF8 80A99DA8 44093000 */  mfc1    $t1, $f6                   
/* 02FFC 80A99DAC 44813000 */  mtc1    $at, $f6                   ## $f6 = 2147483648.00
/* 03000 80A99DB0 24090001 */  addiu   $t1, $zero, 0x0001         ## $t1 = 00000001
/* 03004 80A99DB4 46062181 */  sub.s   $f6, $f4, $f6              
/* 03008 80A99DB8 44C9F800 */  ctc1    $t1, $31
/* 0300C 80A99DBC 00000000 */  nop
/* 03010 80A99DC0 460031A4 */  cvt.w.s $f6, $f6                   
/* 03014 80A99DC4 4449F800 */  cfc1    $t1, $31
/* 03018 80A99DC8 00000000 */  nop
/* 0301C 80A99DCC 31290078 */  andi    $t1, $t1, 0x0078           ## $t1 = 00000000
/* 03020 80A99DD0 15200005 */  bne     $t1, $zero, .L80A99DE8     
/* 03024 80A99DD4 00000000 */  nop
/* 03028 80A99DD8 44093000 */  mfc1    $t1, $f6                   
/* 0302C 80A99DDC 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 03030 80A99DE0 10000007 */  beq     $zero, $zero, .L80A99E00   
/* 03034 80A99DE4 01214825 */  or      $t1, $t1, $at              ## $t1 = 80000000
.L80A99DE8:
/* 03038 80A99DE8 10000005 */  beq     $zero, $zero, .L80A99E00   
/* 0303C 80A99DEC 2409FFFF */  addiu   $t1, $zero, 0xFFFF         ## $t1 = FFFFFFFF
/* 03040 80A99DF0 44093000 */  mfc1    $t1, $f6                   
.L80A99DF4:
/* 03044 80A99DF4 00000000 */  nop
/* 03048 80A99DF8 0520FFFB */  bltz    $t1, .L80A99DE8            
/* 0304C 80A99DFC 00000000 */  nop
.L80A99E00:
/* 03050 80A99E00 A20900C8 */  sb      $t1, 0x00C8($s0)           ## 000000C8
/* 03054 80A99E04 8E450000 */  lw      $a1, 0x0000($s2)           ## 00000000
/* 03058 80A99E08 44C8F800 */  ctc1    $t0, $31
/* 0305C 80A99E0C 0C031AB1 */  jal     Graph_OpenDisp              
/* 03060 80A99E10 00A08825 */  or      $s1, $a1, $zero            ## $s1 = 00000000
/* 03064 80A99E14 C6000220 */  lwc1    $f0, 0x0220($s0)           ## 00000220
/* 03068 80A99E18 240100FF */  addiu   $at, $zero, 0x00FF         ## $at = 000000FF
/* 0306C 80A99E1C 3C0DDB06 */  lui     $t5, 0xDB06                ## $t5 = DB060000
/* 03070 80A99E20 4600020D */  trunc.w.s $f8, $f0                   
/* 03074 80A99E24 44024000 */  mfc1    $v0, $f8                   
/* 03078 80A99E28 00000000 */  nop
/* 0307C 80A99E2C 00021400 */  sll     $v0, $v0, 16               
/* 03080 80A99E30 00021403 */  sra     $v0, $v0, 16               
/* 03084 80A99E34 1441002E */  bne     $v0, $at, .L80A99EF0       
/* 03088 80A99E38 00000000 */  nop
/* 0308C 80A99E3C 8E2202C0 */  lw      $v0, 0x02C0($s1)           ## 000002C0
/* 03090 80A99E40 35AD0020 */  ori     $t5, $t5, 0x0020           ## $t5 = DB060020
/* 03094 80A99E44 240F00FF */  addiu   $t7, $zero, 0x00FF         ## $t7 = 000000FF
/* 03098 80A99E48 244C0008 */  addiu   $t4, $v0, 0x0008           ## $t4 = 00000008
/* 0309C 80A99E4C AE2C02C0 */  sw      $t4, 0x02C0($s1)           ## 000002C0
/* 030A0 80A99E50 AC4D0000 */  sw      $t5, 0x0000($v0)           ## 00000000
/* 030A4 80A99E54 8E440000 */  lw      $a0, 0x0000($s2)           ## 00000000
/* 030A8 80A99E58 AFAF0010 */  sw      $t7, 0x0010($sp)           
/* 030AC 80A99E5C 93A70062 */  lbu     $a3, 0x0062($sp)           
/* 030B0 80A99E60 93A60061 */  lbu     $a2, 0x0061($sp)           
/* 030B4 80A99E64 93A50060 */  lbu     $a1, 0x0060($sp)           
/* 030B8 80A99E68 0C2A6725 */  jal     func_80A99C94              
/* 030BC 80A99E6C AFA20044 */  sw      $v0, 0x0044($sp)           
/* 030C0 80A99E70 8FA30044 */  lw      $v1, 0x0044($sp)           
/* 030C4 80A99E74 3C18DB06 */  lui     $t8, 0xDB06                ## $t8 = DB060000
/* 030C8 80A99E78 37180024 */  ori     $t8, $t8, 0x0024           ## $t8 = DB060024
/* 030CC 80A99E7C AC620004 */  sw      $v0, 0x0004($v1)           ## 00000004
/* 030D0 80A99E80 8E2202C0 */  lw      $v0, 0x02C0($s1)           ## 000002C0
/* 030D4 80A99E84 241900FF */  addiu   $t9, $zero, 0x00FF         ## $t9 = 000000FF
/* 030D8 80A99E88 244B0008 */  addiu   $t3, $v0, 0x0008           ## $t3 = 00000008
/* 030DC 80A99E8C AE2B02C0 */  sw      $t3, 0x02C0($s1)           ## 000002C0
/* 030E0 80A99E90 AC580000 */  sw      $t8, 0x0000($v0)           ## 00000000
/* 030E4 80A99E94 8E440000 */  lw      $a0, 0x0000($s2)           ## 00000000
/* 030E8 80A99E98 AFB90010 */  sw      $t9, 0x0010($sp)           
/* 030EC 80A99E9C 93A7005E */  lbu     $a3, 0x005E($sp)           
/* 030F0 80A99EA0 93A6005D */  lbu     $a2, 0x005D($sp)           
/* 030F4 80A99EA4 93A5005C */  lbu     $a1, 0x005C($sp)           
/* 030F8 80A99EA8 0C2A6725 */  jal     func_80A99C94              
/* 030FC 80A99EAC AFA20040 */  sw      $v0, 0x0040($sp)           
/* 03100 80A99EB0 8FA30040 */  lw      $v1, 0x0040($sp)           
/* 03104 80A99EB4 3C0680AA */  lui     $a2, %hi(func_80A99864)    ## $a2 = 80AA0000
/* 03108 80A99EB8 3C0780AA */  lui     $a3, %hi(func_80A99BC4)    ## $a3 = 80AA0000
/* 0310C 80A99EBC AC620004 */  sw      $v0, 0x0004($v1)           ## 00000004
/* 03110 80A99EC0 AFB00010 */  sw      $s0, 0x0010($sp)           
/* 03114 80A99EC4 C60A0220 */  lwc1    $f10, 0x0220($s0)          ## 00000220
/* 03118 80A99EC8 24E79BC4 */  addiu   $a3, $a3, %lo(func_80A99BC4) ## $a3 = 80A99BC4
/* 0311C 80A99ECC 24C69864 */  addiu   $a2, $a2, %lo(func_80A99864) ## $a2 = 80A99864
/* 03120 80A99ED0 4600540D */  trunc.w.s $f16, $f10                 
/* 03124 80A99ED4 02402025 */  or      $a0, $s2, $zero            ## $a0 = 00000000
/* 03128 80A99ED8 2605014C */  addiu   $a1, $s0, 0x014C           ## $a1 = 0000014C
/* 0312C 80A99EDC 44088000 */  mfc1    $t0, $f16                  
/* 03130 80A99EE0 0C00D2E8 */  jal     func_80034BA0              
/* 03134 80A99EE4 AFA80014 */  sw      $t0, 0x0014($sp)           
/* 03138 80A99EE8 10000071 */  beq     $zero, $zero, .L80A9A0B0   
/* 0313C 80A99EEC 00000000 */  nop
.L80A99EF0:
/* 03140 80A99EF0 1040006F */  beq     $v0, $zero, .L80A9A0B0     
/* 03144 80A99EF4 240A0001 */  addiu   $t2, $zero, 0x0001         ## $t2 = 00000001
/* 03148 80A99EF8 4449F800 */  cfc1    $t1, $31
/* 0314C 80A99EFC 44CAF800 */  ctc1    $t2, $31
/* 03150 80A99F00 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 03154 80A99F04 3C0BDB06 */  lui     $t3, 0xDB06                ## $t3 = DB060000
/* 03158 80A99F08 460004A4 */  cvt.w.s $f18, $f0                  
/* 0315C 80A99F0C 444AF800 */  cfc1    $t2, $31
/* 03160 80A99F10 00000000 */  nop
/* 03164 80A99F14 314A0078 */  andi    $t2, $t2, 0x0078           ## $t2 = 00000000
/* 03168 80A99F18 51400013 */  beql    $t2, $zero, .L80A99F68     
/* 0316C 80A99F1C 440A9000 */  mfc1    $t2, $f18                  
/* 03170 80A99F20 44819000 */  mtc1    $at, $f18                  ## $f18 = 2147483648.00
/* 03174 80A99F24 240A0001 */  addiu   $t2, $zero, 0x0001         ## $t2 = 00000001
/* 03178 80A99F28 46120481 */  sub.s   $f18, $f0, $f18            
/* 0317C 80A99F2C 44CAF800 */  ctc1    $t2, $31
/* 03180 80A99F30 00000000 */  nop
/* 03184 80A99F34 460094A4 */  cvt.w.s $f18, $f18                 
/* 03188 80A99F38 444AF800 */  cfc1    $t2, $31
/* 0318C 80A99F3C 00000000 */  nop
/* 03190 80A99F40 314A0078 */  andi    $t2, $t2, 0x0078           ## $t2 = 00000000
/* 03194 80A99F44 15400005 */  bne     $t2, $zero, .L80A99F5C     
/* 03198 80A99F48 00000000 */  nop
/* 0319C 80A99F4C 440A9000 */  mfc1    $t2, $f18                  
/* 031A0 80A99F50 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 031A4 80A99F54 10000007 */  beq     $zero, $zero, .L80A99F74   
/* 031A8 80A99F58 01415025 */  or      $t2, $t2, $at              ## $t2 = 80000000
.L80A99F5C:
/* 031AC 80A99F5C 10000005 */  beq     $zero, $zero, .L80A99F74   
/* 031B0 80A99F60 240AFFFF */  addiu   $t2, $zero, 0xFFFF         ## $t2 = FFFFFFFF
/* 031B4 80A99F64 440A9000 */  mfc1    $t2, $f18                  
.L80A99F68:
/* 031B8 80A99F68 00000000 */  nop
/* 031BC 80A99F6C 0540FFFB */  bltz    $t2, .L80A99F5C            
/* 031C0 80A99F70 00000000 */  nop
.L80A99F74:
/* 031C4 80A99F74 44C9F800 */  ctc1    $t1, $31
/* 031C8 80A99F78 A3AA0063 */  sb      $t2, 0x0063($sp)           
/* 031CC 80A99F7C 240D0001 */  addiu   $t5, $zero, 0x0001         ## $t5 = 00000001
/* 031D0 80A99F80 C6040220 */  lwc1    $f4, 0x0220($s0)           ## 00000220
/* 031D4 80A99F84 444CF800 */  cfc1    $t4, $31
/* 031D8 80A99F88 44CDF800 */  ctc1    $t5, $31
/* 031DC 80A99F8C 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 031E0 80A99F90 460021A4 */  cvt.w.s $f6, $f4                   
/* 031E4 80A99F94 444DF800 */  cfc1    $t5, $31
/* 031E8 80A99F98 00000000 */  nop
/* 031EC 80A99F9C 31AD0078 */  andi    $t5, $t5, 0x0078           ## $t5 = 00000000
/* 031F0 80A99FA0 51A00013 */  beql    $t5, $zero, .L80A99FF0     
/* 031F4 80A99FA4 440D3000 */  mfc1    $t5, $f6                   
/* 031F8 80A99FA8 44813000 */  mtc1    $at, $f6                   ## $f6 = 2147483648.00
/* 031FC 80A99FAC 240D0001 */  addiu   $t5, $zero, 0x0001         ## $t5 = 00000001
/* 03200 80A99FB0 46062181 */  sub.s   $f6, $f4, $f6              
/* 03204 80A99FB4 44CDF800 */  ctc1    $t5, $31
/* 03208 80A99FB8 00000000 */  nop
/* 0320C 80A99FBC 460031A4 */  cvt.w.s $f6, $f6                   
/* 03210 80A99FC0 444DF800 */  cfc1    $t5, $31
/* 03214 80A99FC4 00000000 */  nop
/* 03218 80A99FC8 31AD0078 */  andi    $t5, $t5, 0x0078           ## $t5 = 00000000
/* 0321C 80A99FCC 15A00005 */  bne     $t5, $zero, .L80A99FE4     
/* 03220 80A99FD0 00000000 */  nop
/* 03224 80A99FD4 440D3000 */  mfc1    $t5, $f6                   
/* 03228 80A99FD8 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 0322C 80A99FDC 10000007 */  beq     $zero, $zero, .L80A99FFC   
/* 03230 80A99FE0 01A16825 */  or      $t5, $t5, $at              ## $t5 = 80000000
.L80A99FE4:
/* 03234 80A99FE4 10000005 */  beq     $zero, $zero, .L80A99FFC   
/* 03238 80A99FE8 240DFFFF */  addiu   $t5, $zero, 0xFFFF         ## $t5 = FFFFFFFF
/* 0323C 80A99FEC 440D3000 */  mfc1    $t5, $f6                   
.L80A99FF0:
/* 03240 80A99FF0 00000000 */  nop
/* 03244 80A99FF4 05A0FFFB */  bltz    $t5, .L80A99FE4            
/* 03248 80A99FF8 00000000 */  nop
.L80A99FFC:
/* 0324C 80A99FFC A3AD005F */  sb      $t5, 0x005F($sp)           
/* 03250 80A9A000 8E2202D0 */  lw      $v0, 0x02D0($s1)           ## 000002D0
/* 03254 80A9A004 356B0020 */  ori     $t3, $t3, 0x0020           ## $t3 = DB060020
/* 03258 80A9A008 44CCF800 */  ctc1    $t4, $31
/* 0325C 80A9A00C 244F0008 */  addiu   $t7, $v0, 0x0008           ## $t7 = 00000008
/* 03260 80A9A010 AE2F02D0 */  sw      $t7, 0x02D0($s1)           ## 000002D0
/* 03264 80A9A014 AC4B0000 */  sw      $t3, 0x0000($v0)           ## 00000000
/* 03268 80A9A018 93B80063 */  lbu     $t8, 0x0063($sp)           
/* 0326C 80A9A01C 8E440000 */  lw      $a0, 0x0000($s2)           ## 00000000
/* 03270 80A9A020 93A70062 */  lbu     $a3, 0x0062($sp)           
/* 03274 80A9A024 93A60061 */  lbu     $a2, 0x0061($sp)           
/* 03278 80A9A028 93A50060 */  lbu     $a1, 0x0060($sp)           
/* 0327C 80A9A02C AFA2003C */  sw      $v0, 0x003C($sp)           
/* 03280 80A9A030 0C2A6725 */  jal     func_80A99C94              
/* 03284 80A9A034 AFB80010 */  sw      $t8, 0x0010($sp)           
/* 03288 80A9A038 8FA3003C */  lw      $v1, 0x003C($sp)           
/* 0328C 80A9A03C 3C0EDB06 */  lui     $t6, 0xDB06                ## $t6 = DB060000
/* 03290 80A9A040 35CE0024 */  ori     $t6, $t6, 0x0024           ## $t6 = DB060024
/* 03294 80A9A044 AC620004 */  sw      $v0, 0x0004($v1)           ## 00000004
/* 03298 80A9A048 8E2202D0 */  lw      $v0, 0x02D0($s1)           ## 000002D0
/* 0329C 80A9A04C 24590008 */  addiu   $t9, $v0, 0x0008           ## $t9 = 00000008
/* 032A0 80A9A050 AE3902D0 */  sw      $t9, 0x02D0($s1)           ## 000002D0
/* 032A4 80A9A054 AC4E0000 */  sw      $t6, 0x0000($v0)           ## 00000000
/* 032A8 80A9A058 93A8005F */  lbu     $t0, 0x005F($sp)           
/* 032AC 80A9A05C 8E440000 */  lw      $a0, 0x0000($s2)           ## 00000000
/* 032B0 80A9A060 93A7005E */  lbu     $a3, 0x005E($sp)           
/* 032B4 80A9A064 93A6005D */  lbu     $a2, 0x005D($sp)           
/* 032B8 80A9A068 93A5005C */  lbu     $a1, 0x005C($sp)           
/* 032BC 80A9A06C AFA20038 */  sw      $v0, 0x0038($sp)           
/* 032C0 80A9A070 0C2A6725 */  jal     func_80A99C94              
/* 032C4 80A9A074 AFA80010 */  sw      $t0, 0x0010($sp)           
/* 032C8 80A9A078 8FA30038 */  lw      $v1, 0x0038($sp)           
/* 032CC 80A9A07C 3C0680AA */  lui     $a2, %hi(func_80A99864)    ## $a2 = 80AA0000
/* 032D0 80A9A080 3C0780AA */  lui     $a3, %hi(func_80A99BC4)    ## $a3 = 80AA0000
/* 032D4 80A9A084 AC620004 */  sw      $v0, 0x0004($v1)           ## 00000004
/* 032D8 80A9A088 AFB00010 */  sw      $s0, 0x0010($sp)           
/* 032DC 80A9A08C C6080220 */  lwc1    $f8, 0x0220($s0)           ## 00000220
/* 032E0 80A9A090 24E79BC4 */  addiu   $a3, $a3, %lo(func_80A99BC4) ## $a3 = 80A99BC4
/* 032E4 80A9A094 24C69864 */  addiu   $a2, $a2, %lo(func_80A99864) ## $a2 = 80A99864
/* 032E8 80A9A098 4600428D */  trunc.w.s $f10, $f8                  
/* 032EC 80A9A09C 02402025 */  or      $a0, $s2, $zero            ## $a0 = 00000000
/* 032F0 80A9A0A0 2605014C */  addiu   $a1, $s0, 0x014C           ## $a1 = 0000014C
/* 032F4 80A9A0A4 440A5000 */  mfc1    $t2, $f10                  
/* 032F8 80A9A0A8 0C00D331 */  jal     func_80034CC4              
/* 032FC 80A9A0AC AFAA0014 */  sw      $t2, 0x0014($sp)           
.L80A9A0B0:
/* 03300 80A9A0B0 3C0680AA */  lui     $a2, %hi(D_80A9A7AC)       ## $a2 = 80AA0000
/* 03304 80A9A0B4 24C6A7AC */  addiu   $a2, $a2, %lo(D_80A9A7AC)  ## $a2 = 80A9A7AC
/* 03308 80A9A0B8 27A40048 */  addiu   $a0, $sp, 0x0048           ## $a0 = FFFFFFE0
/* 0330C 80A9A0BC 8E450000 */  lw      $a1, 0x0000($s2)           ## 00000000
/* 03310 80A9A0C0 0C031AD5 */  jal     Graph_CloseDisp              
/* 03314 80A9A0C4 24070858 */  addiu   $a3, $zero, 0x0858         ## $a3 = 00000858
/* 03318 80A9A0C8 8FBF002C */  lw      $ra, 0x002C($sp)           
/* 0331C 80A9A0CC 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 03320 80A9A0D0 8FB10024 */  lw      $s1, 0x0024($sp)           
/* 03324 80A9A0D4 8FB20028 */  lw      $s2, 0x0028($sp)           
/* 03328 80A9A0D8 03E00008 */  jr      $ra                        
/* 0332C 80A9A0DC 27BD0068 */  addiu   $sp, $sp, 0x0068           ## $sp = 00000000

