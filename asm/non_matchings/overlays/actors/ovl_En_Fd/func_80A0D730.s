glabel func_80A0D730
/* 02070 80A0D730 27BDFF50 */  addiu   $sp, $sp, 0xFF50           ## $sp = FFFFFF50
/* 02074 80A0D734 AFBF005C */  sw      $ra, 0x005C($sp)           
/* 02078 80A0D738 AFBE0058 */  sw      $s8, 0x0058($sp)           
/* 0207C 80A0D73C AFB70054 */  sw      $s7, 0x0054($sp)           
/* 02080 80A0D740 AFB60050 */  sw      $s6, 0x0050($sp)           
/* 02084 80A0D744 AFB5004C */  sw      $s5, 0x004C($sp)           
/* 02088 80A0D748 AFB40048 */  sw      $s4, 0x0048($sp)           
/* 0208C 80A0D74C AFB30044 */  sw      $s3, 0x0044($sp)           
/* 02090 80A0D750 AFB20040 */  sw      $s2, 0x0040($sp)           
/* 02094 80A0D754 AFB1003C */  sw      $s1, 0x003C($sp)           
/* 02098 80A0D758 AFB00038 */  sw      $s0, 0x0038($sp)           
/* 0209C 80A0D75C F7BA0030 */  sdc1    $f26, 0x0030($sp)          
/* 020A0 80A0D760 F7B80028 */  sdc1    $f24, 0x0028($sp)          
/* 020A4 80A0D764 F7B60020 */  sdc1    $f22, 0x0020($sp)          
/* 020A8 80A0D768 F7B40018 */  sdc1    $f20, 0x0018($sp)          
/* 020AC 80A0D76C AFA500B4 */  sw      $a1, 0x00B4($sp)           
/* 020B0 80A0D770 8CA50000 */  lw      $a1, 0x0000($a1)           ## 00000000
/* 020B4 80A0D774 0080B825 */  or      $s7, $a0, $zero            ## $s7 = 00000000
/* 020B8 80A0D778 24920620 */  addiu   $s2, $a0, 0x0620           ## $s2 = 00000620
/* 020BC 80A0D77C 3C0680A1 */  lui     $a2, %hi(D_80A0E140)       ## $a2 = 80A10000
/* 020C0 80A0D780 24C6E140 */  addiu   $a2, $a2, %lo(D_80A0E140)  ## $a2 = 80A0E140
/* 020C4 80A0D784 27A40090 */  addiu   $a0, $sp, 0x0090           ## $a0 = FFFFFFE0
/* 020C8 80A0D788 240707B1 */  addiu   $a3, $zero, 0x07B1         ## $a3 = 000007B1
/* 020CC 80A0D78C 0C031AB1 */  jal     func_800C6AC4              
/* 020D0 80A0D790 00A08025 */  or      $s0, $a1, $zero            ## $s0 = 00000000
/* 020D4 80A0D794 8FAF00B4 */  lw      $t7, 0x00B4($sp)           
/* 020D8 80A0D798 0000F025 */  or      $s8, $zero, $zero          ## $s8 = 00000000
/* 020DC 80A0D79C 0C024F61 */  jal     func_80093D84              
/* 020E0 80A0D7A0 8DE40000 */  lw      $a0, 0x0000($t7)           ## 00000000
/* 020E4 80A0D7A4 3C014100 */  lui     $at, 0x4100                ## $at = 41000000
/* 020E8 80A0D7A8 4481D000 */  mtc1    $at, $f26                  ## $f26 = 8.00
/* 020EC 80A0D7AC 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 020F0 80A0D7B0 4481C000 */  mtc1    $at, $f24                  ## $f24 = 1.00
/* 020F4 80A0D7B4 3C01437F */  lui     $at, 0x437F                ## $at = 437F0000
/* 020F8 80A0D7B8 4481B000 */  mtc1    $at, $f22                  ## $f22 = 255.00
/* 020FC 80A0D7BC 4481A000 */  mtc1    $at, $f20                  ## $f20 = 255.00
/* 02100 80A0D7C0 0000A025 */  or      $s4, $zero, $zero          ## $s4 = 00000000
/* 02104 80A0D7C4 3C16DE00 */  lui     $s6, 0xDE00                ## $s6 = DE000000
.L80A0D7C8:
/* 02108 80A0D7C8 92580000 */  lbu     $t8, 0x0000($s2)           ## 00000620
/* 0210C 80A0D7CC 24010001 */  addiu   $at, $zero, 0x0001         ## $at = 00000001
/* 02110 80A0D7D0 8FB300B4 */  lw      $s3, 0x00B4($sp)           
/* 02114 80A0D7D4 170100BE */  bne     $t8, $at, .L80A0DAD0       
/* 02118 80A0D7D8 3C150600 */  lui     $s5, 0x0600                ## $s5 = 06000000
/* 0211C 80A0D7DC 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 02120 80A0D7E0 34211DA0 */  ori     $at, $at, 0x1DA0           ## $at = 00011DA0
/* 02124 80A0D7E4 02619821 */  addu    $s3, $s3, $at              
/* 02128 80A0D7E8 17C00039 */  bne     $s8, $zero, .L80A0D8D0     
/* 0212C 80A0D7EC 26B57938 */  addiu   $s5, $s5, 0x7938           ## $s5 = 06007938
/* 02130 80A0D7F0 3C110600 */  lui     $s1, 0x0600                ## $s1 = 06000000
/* 02134 80A0D7F4 26317928 */  addiu   $s1, $s1, 0x7928           ## $s1 = 06007928
/* 02138 80A0D7F8 8E0402D0 */  lw      $a0, 0x02D0($s0)           ## 000002D0
/* 0213C 80A0D7FC 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 02140 80A0D800 0C024DDD */  jal     Gfx_CallSetupDL              
/* 02144 80A0D804 241E0001 */  addiu   $s8, $zero, 0x0001         ## $s8 = 00000001
/* 02148 80A0D808 AE0202D0 */  sw      $v0, 0x02D0($s0)           ## 000002D0
/* 0214C 80A0D80C 24590008 */  addiu   $t9, $v0, 0x0008           ## $t9 = 00000008
/* 02150 80A0D810 AE1902D0 */  sw      $t9, 0x02D0($s0)           ## 000002D0
/* 02154 80A0D814 AC510004 */  sw      $s1, 0x0004($v0)           ## 00000004
/* 02158 80A0D818 AC560000 */  sw      $s6, 0x0000($v0)           ## 00000000
/* 0215C 80A0D81C 8E0302D0 */  lw      $v1, 0x02D0($s0)           ## 000002D0
/* 02160 80A0D820 3C09FB00 */  lui     $t1, 0xFB00                ## $t1 = FB000000
/* 02164 80A0D824 240B0001 */  addiu   $t3, $zero, 0x0001         ## $t3 = 00000001
/* 02168 80A0D828 24680008 */  addiu   $t0, $v1, 0x0008           ## $t0 = 00000008
/* 0216C 80A0D82C AE0802D0 */  sw      $t0, 0x02D0($s0)           ## 000002D0
/* 02170 80A0D830 AC690000 */  sw      $t1, 0x0000($v1)           ## 00000000
/* 02174 80A0D834 C6E404CC */  lwc1    $f4, 0x04CC($s7)           ## 000004CC
/* 02178 80A0D838 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 0217C 80A0D83C 46142183 */  div.s   $f6, $f4, $f20             
/* 02180 80A0D840 46163202 */  mul.s   $f8, $f6, $f22             
/* 02184 80A0D844 444AF800 */  cfc1    $t2, $31
/* 02188 80A0D848 44CBF800 */  ctc1    $t3, $31
/* 0218C 80A0D84C 00000000 */  nop
/* 02190 80A0D850 460042A4 */  cvt.w.s $f10, $f8                  
/* 02194 80A0D854 444BF800 */  cfc1    $t3, $31
/* 02198 80A0D858 00000000 */  nop
/* 0219C 80A0D85C 316B0078 */  andi    $t3, $t3, 0x0078           ## $t3 = 00000000
/* 021A0 80A0D860 51600013 */  beql    $t3, $zero, .L80A0D8B0     
/* 021A4 80A0D864 440B5000 */  mfc1    $t3, $f10                  
/* 021A8 80A0D868 44815000 */  mtc1    $at, $f10                  ## $f10 = 2147483648.00
/* 021AC 80A0D86C 240B0001 */  addiu   $t3, $zero, 0x0001         ## $t3 = 00000001
/* 021B0 80A0D870 460A4281 */  sub.s   $f10, $f8, $f10            
/* 021B4 80A0D874 44CBF800 */  ctc1    $t3, $31
/* 021B8 80A0D878 00000000 */  nop
/* 021BC 80A0D87C 460052A4 */  cvt.w.s $f10, $f10                 
/* 021C0 80A0D880 444BF800 */  cfc1    $t3, $31
/* 021C4 80A0D884 00000000 */  nop
/* 021C8 80A0D888 316B0078 */  andi    $t3, $t3, 0x0078           ## $t3 = 00000000
/* 021CC 80A0D88C 15600005 */  bne     $t3, $zero, .L80A0D8A4     
/* 021D0 80A0D890 00000000 */  nop
/* 021D4 80A0D894 440B5000 */  mfc1    $t3, $f10                  
/* 021D8 80A0D898 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 021DC 80A0D89C 10000007 */  beq     $zero, $zero, .L80A0D8BC   
/* 021E0 80A0D8A0 01615825 */  or      $t3, $t3, $at              ## $t3 = 80000000
.L80A0D8A4:
/* 021E4 80A0D8A4 10000005 */  beq     $zero, $zero, .L80A0D8BC   
/* 021E8 80A0D8A8 240BFFFF */  addiu   $t3, $zero, 0xFFFF         ## $t3 = FFFFFFFF
/* 021EC 80A0D8AC 440B5000 */  mfc1    $t3, $f10                  
.L80A0D8B0:
/* 021F0 80A0D8B0 00000000 */  nop
/* 021F4 80A0D8B4 0560FFFB */  bltz    $t3, .L80A0D8A4            
/* 021F8 80A0D8B8 00000000 */  nop
.L80A0D8BC:
/* 021FC 80A0D8BC 316D00FF */  andi    $t5, $t3, 0x00FF           ## $t5 = 000000FF
/* 02200 80A0D8C0 3C01FF0A */  lui     $at, 0xFF0A                ## $at = FF0A0000
/* 02204 80A0D8C4 44CAF800 */  ctc1    $t2, $31
/* 02208 80A0D8C8 01A17025 */  or      $t6, $t5, $at              ## $t6 = FF0A00FF
/* 0220C 80A0D8CC AC6E0004 */  sw      $t6, 0x0004($v1)           ## 00000004
.L80A0D8D0:
/* 02210 80A0D8D0 8E0302D0 */  lw      $v1, 0x02D0($s0)           ## 000002D0
/* 02214 80A0D8D4 3C18FA00 */  lui     $t8, 0xFA00                ## $t8 = FA000000
/* 02218 80A0D8D8 24080001 */  addiu   $t0, $zero, 0x0001         ## $t0 = 00000001
/* 0221C 80A0D8DC 246F0008 */  addiu   $t7, $v1, 0x0008           ## $t7 = 00000008
/* 02220 80A0D8E0 AE0F02D0 */  sw      $t7, 0x02D0($s0)           ## 000002D0
/* 02224 80A0D8E4 AC780000 */  sw      $t8, 0x0000($v1)           ## 00000000
/* 02228 80A0D8E8 C6F004CC */  lwc1    $f16, 0x04CC($s7)          ## 000004CC
/* 0222C 80A0D8EC 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 02230 80A0D8F0 3C0DE700 */  lui     $t5, 0xE700                ## $t5 = E7000000
/* 02234 80A0D8F4 46148483 */  div.s   $f18, $f16, $f20           
/* 02238 80A0D8F8 46169102 */  mul.s   $f4, $f18, $f22            
/* 0223C 80A0D8FC 4459F800 */  cfc1    $t9, $31
/* 02240 80A0D900 44C8F800 */  ctc1    $t0, $31
/* 02244 80A0D904 00000000 */  nop
/* 02248 80A0D908 460021A4 */  cvt.w.s $f6, $f4                   
/* 0224C 80A0D90C 4448F800 */  cfc1    $t0, $31
/* 02250 80A0D910 00000000 */  nop
/* 02254 80A0D914 31080078 */  andi    $t0, $t0, 0x0078           ## $t0 = 00000000
/* 02258 80A0D918 51000013 */  beql    $t0, $zero, .L80A0D968     
/* 0225C 80A0D91C 44083000 */  mfc1    $t0, $f6                   
/* 02260 80A0D920 44813000 */  mtc1    $at, $f6                   ## $f6 = 2147483648.00
/* 02264 80A0D924 24080001 */  addiu   $t0, $zero, 0x0001         ## $t0 = 00000001
/* 02268 80A0D928 46062181 */  sub.s   $f6, $f4, $f6              
/* 0226C 80A0D92C 44C8F800 */  ctc1    $t0, $31
/* 02270 80A0D930 00000000 */  nop
/* 02274 80A0D934 460031A4 */  cvt.w.s $f6, $f6                   
/* 02278 80A0D938 4448F800 */  cfc1    $t0, $31
/* 0227C 80A0D93C 00000000 */  nop
/* 02280 80A0D940 31080078 */  andi    $t0, $t0, 0x0078           ## $t0 = 00000000
/* 02284 80A0D944 15000005 */  bne     $t0, $zero, .L80A0D95C     
/* 02288 80A0D948 00000000 */  nop
/* 0228C 80A0D94C 44083000 */  mfc1    $t0, $f6                   
/* 02290 80A0D950 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 02294 80A0D954 10000007 */  beq     $zero, $zero, .L80A0D974   
/* 02298 80A0D958 01014025 */  or      $t0, $t0, $at              ## $t0 = 80000000
.L80A0D95C:
/* 0229C 80A0D95C 10000005 */  beq     $zero, $zero, .L80A0D974   
/* 022A0 80A0D960 2408FFFF */  addiu   $t0, $zero, 0xFFFF         ## $t0 = FFFFFFFF
/* 022A4 80A0D964 44083000 */  mfc1    $t0, $f6                   
.L80A0D968:
/* 022A8 80A0D968 00000000 */  nop
/* 022AC 80A0D96C 0500FFFB */  bltz    $t0, .L80A0D95C            
/* 022B0 80A0D970 00000000 */  nop
.L80A0D974:
/* 022B4 80A0D974 310A00FF */  andi    $t2, $t0, 0x00FF           ## $t2 = 000000FF
/* 022B8 80A0D978 3C01FFFF */  lui     $at, 0xFFFF                ## $at = FFFF0000
/* 022BC 80A0D97C 01415825 */  or      $t3, $t2, $at              ## $t3 = FFFF00FF
/* 022C0 80A0D980 AC6B0004 */  sw      $t3, 0x0004($v1)           ## 00000004
/* 022C4 80A0D984 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 022C8 80A0D988 44D9F800 */  ctc1    $t9, $31
/* 022CC 80A0D98C 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 022D0 80A0D990 244C0008 */  addiu   $t4, $v0, 0x0008           ## $t4 = 00000008
/* 022D4 80A0D994 AE0C02D0 */  sw      $t4, 0x02D0($s0)           ## 000002D0
/* 022D8 80A0D998 AC400004 */  sw      $zero, 0x0004($v0)         ## 00000004
/* 022DC 80A0D99C AC4D0000 */  sw      $t5, 0x0000($v0)           ## 00000000
/* 022E0 80A0D9A0 8E46001C */  lw      $a2, 0x001C($s2)           ## 0000063C
/* 022E4 80A0D9A4 C64E0018 */  lwc1    $f14, 0x0018($s2)          ## 00000638
/* 022E8 80A0D9A8 0C034261 */  jal     Matrix_Translate              
/* 022EC 80A0D9AC C64C0014 */  lwc1    $f12, 0x0014($s2)          ## 00000634
/* 022F0 80A0D9B0 0C0347F5 */  jal     func_800D1FD4              
/* 022F4 80A0D9B4 02602025 */  or      $a0, $s3, $zero            ## $a0 = 00000000
/* 022F8 80A0D9B8 C64C0004 */  lwc1    $f12, 0x0004($s2)          ## 00000624
/* 022FC 80A0D9BC 4406C000 */  mfc1    $a2, $f24                  
/* 02300 80A0D9C0 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 02304 80A0D9C4 0C0342A3 */  jal     Matrix_Scale              
/* 02308 80A0D9C8 46006386 */  mov.s   $f14, $f12                 
/* 0230C 80A0D9CC 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 02310 80A0D9D0 3C0FDA38 */  lui     $t7, 0xDA38                ## $t7 = DA380000
/* 02314 80A0D9D4 35EF0003 */  ori     $t7, $t7, 0x0003           ## $t7 = DA380003
/* 02318 80A0D9D8 244E0008 */  addiu   $t6, $v0, 0x0008           ## $t6 = 00000008
/* 0231C 80A0D9DC AE0E02D0 */  sw      $t6, 0x02D0($s0)           ## 000002D0
/* 02320 80A0D9E0 AC4F0000 */  sw      $t7, 0x0000($v0)           ## 00000000
/* 02324 80A0D9E4 8FB800B4 */  lw      $t8, 0x00B4($sp)           
/* 02328 80A0D9E8 3C0580A1 */  lui     $a1, %hi(D_80A0E150)       ## $a1 = 80A10000
/* 0232C 80A0D9EC 24A5E150 */  addiu   $a1, $a1, %lo(D_80A0E150)  ## $a1 = 80A0E150
/* 02330 80A0D9F0 240607D6 */  addiu   $a2, $zero, 0x07D6         ## $a2 = 000007D6
/* 02334 80A0D9F4 00408825 */  or      $s1, $v0, $zero            ## $s1 = 00000000
/* 02338 80A0D9F8 0C0346A2 */  jal     Matrix_NewMtx              
/* 0233C 80A0D9FC 8F040000 */  lw      $a0, 0x0000($t8)           ## 00000000
/* 02340 80A0DA00 AE220004 */  sw      $v0, 0x0004($s1)           ## 00000004
/* 02344 80A0DA04 92590001 */  lbu     $t9, 0x0001($s2)           ## 00000621
/* 02348 80A0DA08 3C014F80 */  lui     $at, 0x4F80                ## $at = 4F800000
/* 0234C 80A0DA0C 44994000 */  mtc1    $t9, $f8                   ## $f8 = 0.00
/* 02350 80A0DA10 07210004 */  bgez    $t9, .L80A0DA24            
/* 02354 80A0DA14 468042A0 */  cvt.s.w $f10, $f8                  
/* 02358 80A0DA18 44818000 */  mtc1    $at, $f16                  ## $f16 = 4294967296.00
/* 0235C 80A0DA1C 00000000 */  nop
/* 02360 80A0DA20 46105280 */  add.s   $f10, $f10, $f16           
.L80A0DA24:
/* 02364 80A0DA24 92480002 */  lbu     $t0, 0x0002($s2)           ## 00000622
/* 02368 80A0DA28 3C014F80 */  lui     $at, 0x4F80                ## $at = 4F800000
/* 0236C 80A0DA2C 44889000 */  mtc1    $t0, $f18                  ## $f18 = 0.00
/* 02370 80A0DA30 05010004 */  bgez    $t0, .L80A0DA44            
/* 02374 80A0DA34 46809120 */  cvt.s.w $f4, $f18                  
/* 02378 80A0DA38 44813000 */  mtc1    $at, $f6                   ## $f6 = 4294967296.00
/* 0237C 80A0DA3C 00000000 */  nop
/* 02380 80A0DA40 46062100 */  add.s   $f4, $f4, $f6              
.L80A0DA44:
/* 02384 80A0DA44 4604D203 */  div.s   $f8, $f26, $f4             
/* 02388 80A0DA48 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 0238C 80A0DA4C 3C0BDB06 */  lui     $t3, 0xDB06                ## $t3 = DB060000
/* 02390 80A0DA50 356B0020 */  ori     $t3, $t3, 0x0020           ## $t3 = DB060020
/* 02394 80A0DA54 244A0008 */  addiu   $t2, $v0, 0x0008           ## $t2 = 00000008
/* 02398 80A0DA58 AE0A02D0 */  sw      $t2, 0x02D0($s0)           ## 000002D0
/* 0239C 80A0DA5C 3C0480A1 */  lui     $a0, %hi(D_80A0E0F8)       ## $a0 = 80A10000
/* 023A0 80A0DA60 AC4B0000 */  sw      $t3, 0x0000($v0)           ## 00000000
/* 023A4 80A0DA64 3C198016 */  lui     $t9, 0x8016                ## $t9 = 80160000
/* 023A8 80A0DA68 3C0100FF */  lui     $at, 0x00FF                ## $at = 00FF0000
/* 023AC 80A0DA6C 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = 00FFFFFF
/* 023B0 80A0DA70 46085402 */  mul.s   $f16, $f10, $f8            
/* 023B4 80A0DA74 4600848D */  trunc.w.s $f18, $f16                 
/* 023B8 80A0DA78 44059000 */  mfc1    $a1, $f18                  
/* 023BC 80A0DA7C 00000000 */  nop
/* 023C0 80A0DA80 00052C00 */  sll     $a1, $a1, 16               
/* 023C4 80A0DA84 00052C03 */  sra     $a1, $a1, 16               
/* 023C8 80A0DA88 00056080 */  sll     $t4, $a1,  2               
/* 023CC 80A0DA8C 008C2021 */  addu    $a0, $a0, $t4              
/* 023D0 80A0DA90 8C84E0F8 */  lw      $a0, %lo(D_80A0E0F8)($a0)  
/* 023D4 80A0DA94 00047100 */  sll     $t6, $a0,  4               
/* 023D8 80A0DA98 000E7F02 */  srl     $t7, $t6, 28               
/* 023DC 80A0DA9C 000FC080 */  sll     $t8, $t7,  2               
/* 023E0 80A0DAA0 0338C821 */  addu    $t9, $t9, $t8              
/* 023E4 80A0DAA4 8F396FA8 */  lw      $t9, 0x6FA8($t9)           ## 80166FA8
/* 023E8 80A0DAA8 00816824 */  and     $t5, $a0, $at              
/* 023EC 80A0DAAC 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 023F0 80A0DAB0 01B94021 */  addu    $t0, $t5, $t9              
/* 023F4 80A0DAB4 01014821 */  addu    $t1, $t0, $at              
/* 023F8 80A0DAB8 AC490004 */  sw      $t1, 0x0004($v0)           ## 00000004
/* 023FC 80A0DABC 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 02400 80A0DAC0 244A0008 */  addiu   $t2, $v0, 0x0008           ## $t2 = 00000008
/* 02404 80A0DAC4 AE0A02D0 */  sw      $t2, 0x02D0($s0)           ## 000002D0
/* 02408 80A0DAC8 AC550004 */  sw      $s5, 0x0004($v0)           ## 00000004
/* 0240C 80A0DACC AC560000 */  sw      $s6, 0x0000($v0)           ## 00000000
.L80A0DAD0:
/* 02410 80A0DAD0 26940001 */  addiu   $s4, $s4, 0x0001           ## $s4 = 00000001
/* 02414 80A0DAD4 0014A400 */  sll     $s4, $s4, 16               
/* 02418 80A0DAD8 0014A403 */  sra     $s4, $s4, 16               
/* 0241C 80A0DADC 2A8100C8 */  slti    $at, $s4, 0x00C8           
/* 02420 80A0DAE0 1420FF39 */  bne     $at, $zero, .L80A0D7C8     
/* 02424 80A0DAE4 26520038 */  addiu   $s2, $s2, 0x0038           ## $s2 = 00000658
/* 02428 80A0DAE8 8FAB00B4 */  lw      $t3, 0x00B4($sp)           
/* 0242C 80A0DAEC 3C0680A1 */  lui     $a2, %hi(D_80A0E160)       ## $a2 = 80A10000
/* 02430 80A0DAF0 24C6E160 */  addiu   $a2, $a2, %lo(D_80A0E160)  ## $a2 = 80A0E160
/* 02434 80A0DAF4 27A40090 */  addiu   $a0, $sp, 0x0090           ## $a0 = FFFFFFE0
/* 02438 80A0DAF8 240707E4 */  addiu   $a3, $zero, 0x07E4         ## $a3 = 000007E4
/* 0243C 80A0DAFC 0C031AD5 */  jal     func_800C6B54              
/* 02440 80A0DB00 8D650000 */  lw      $a1, 0x0000($t3)           ## DB060020
/* 02444 80A0DB04 8FBF005C */  lw      $ra, 0x005C($sp)           
/* 02448 80A0DB08 D7B40018 */  ldc1    $f20, 0x0018($sp)          
/* 0244C 80A0DB0C D7B60020 */  ldc1    $f22, 0x0020($sp)          
/* 02450 80A0DB10 D7B80028 */  ldc1    $f24, 0x0028($sp)          
/* 02454 80A0DB14 D7BA0030 */  ldc1    $f26, 0x0030($sp)          
/* 02458 80A0DB18 8FB00038 */  lw      $s0, 0x0038($sp)           
/* 0245C 80A0DB1C 8FB1003C */  lw      $s1, 0x003C($sp)           
/* 02460 80A0DB20 8FB20040 */  lw      $s2, 0x0040($sp)           
/* 02464 80A0DB24 8FB30044 */  lw      $s3, 0x0044($sp)           
/* 02468 80A0DB28 8FB40048 */  lw      $s4, 0x0048($sp)           
/* 0246C 80A0DB2C 8FB5004C */  lw      $s5, 0x004C($sp)           
/* 02470 80A0DB30 8FB60050 */  lw      $s6, 0x0050($sp)           
/* 02474 80A0DB34 8FB70054 */  lw      $s7, 0x0054($sp)           
/* 02478 80A0DB38 8FBE0058 */  lw      $s8, 0x0058($sp)           
/* 0247C 80A0DB3C 03E00008 */  jr      $ra                        
/* 02480 80A0DB40 27BD00B0 */  addiu   $sp, $sp, 0x00B0           ## $sp = 00000000


