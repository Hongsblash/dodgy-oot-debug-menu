.late_rodata
glabel D_80867B98
    .float 950
glabel D_80867B9C
    .float 0.33333334
glabel D_80867BA0
    .float 0.041666668
glabel D_80867BA4
    .float 0.1

.text
glabel func_80865F1C
/* 0020C 80865F1C 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 00210 80865F20 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 00214 80865F24 AFA5001C */  sw      $a1, 0x001C($sp)           
/* 00218 80865F28 3C014248 */  lui     $at, 0x4248                ## $at = 42480000
/* 0021C 80865F2C 44811000 */  mtc1    $at, $f2                   ## $f2 = 50.00
/* 00220 80865F30 C48000F0 */  lwc1    $f0, 0x00F0($a0)           ## 000000F0
/* 00224 80865F34 3C018086 */  lui     $at, %hi(D_80867B98)       ## $at = 80860000
/* 00228 80865F38 4602003C */  c.lt.s  $f0, $f2                   
/* 0022C 80865F3C 00000000 */  nop
/* 00230 80865F40 45000005 */  bc1f    .L80865F58                 
/* 00234 80865F44 00000000 */  nop
/* 00238 80865F48 3C014120 */  lui     $at, 0x4120                ## $at = 41200000
/* 0023C 80865F4C 44816000 */  mtc1    $at, $f12                  ## $f12 = 10.00
/* 00240 80865F50 10000012 */  beq     $zero, $zero, .L80865F9C   
/* 00244 80865F54 94830166 */  lhu     $v1, 0x0166($a0)           ## 00000166
.L80865F58:
/* 00248 80865F58 C4247B98 */  lwc1    $f4, %lo(D_80867B98)($at)  
/* 0024C 80865F5C 3C018086 */  lui     $at, %hi(D_80867B9C)       ## $at = 80860000
/* 00250 80865F60 4600203C */  c.lt.s  $f4, $f0                   
/* 00254 80865F64 00000000 */  nop
/* 00258 80865F68 45020006 */  bc1fl   .L80865F84                 
/* 0025C 80865F6C 46020181 */  sub.s   $f6, $f0, $f2              
/* 00260 80865F70 3C01439B */  lui     $at, 0x439B                ## $at = 439B0000
/* 00264 80865F74 44816000 */  mtc1    $at, $f12                  ## $f12 = 310.00
/* 00268 80865F78 10000008 */  beq     $zero, $zero, .L80865F9C   
/* 0026C 80865F7C 94830166 */  lhu     $v1, 0x0166($a0)           ## 00000166
/* 00270 80865F80 46020181 */  sub.s   $f6, $f0, $f2              
.L80865F84:
/* 00274 80865F84 C4287B9C */  lwc1    $f8, %lo(D_80867B9C)($at)  
/* 00278 80865F88 3C014120 */  lui     $at, 0x4120                ## $at = 41200000
/* 0027C 80865F8C 44818000 */  mtc1    $at, $f16                  ## $f16 = 10.00
/* 00280 80865F90 46083282 */  mul.s   $f10, $f6, $f8             
/* 00284 80865F94 46105300 */  add.s   $f12, $f10, $f16           
/* 00288 80865F98 94830166 */  lhu     $v1, 0x0166($a0)           ## 00000166
.L80865F9C:
/* 0028C 80865F9C 240A00FF */  addiu   $t2, $zero, 0x00FF         ## $t2 = 000000FF
/* 00290 80865FA0 10600029 */  beq     $v1, $zero, .L80866048     
/* 00294 80865FA4 246EFFFF */  addiu   $t6, $v1, 0xFFFF           ## $t6 = FFFFFFFF
/* 00298 80865FA8 31C3FFFF */  andi    $v1, $t6, 0xFFFF           ## $v1 = 0000FFFF
/* 0029C 80865FAC 28610008 */  slti    $at, $v1, 0x0008           
/* 002A0 80865FB0 A48E0166 */  sh      $t6, 0x0166($a0)           ## 00000166
/* 002A4 80865FB4 14200024 */  bne     $at, $zero, .L80866048     
/* 002A8 80865FB8 00601025 */  or      $v0, $v1, $zero            ## $v0 = 0000FFFF
/* 002AC 80865FBC 244FFFF8 */  addiu   $t7, $v0, 0xFFF8           ## $t7 = 0000FFF7
/* 002B0 80865FC0 448F9000 */  mtc1    $t7, $f18                  ## $f18 = 0.00
/* 002B4 80865FC4 3C018086 */  lui     $at, %hi(D_80867BA0)       ## $at = 80860000
/* 002B8 80865FC8 C4267BA0 */  lwc1    $f6, %lo(D_80867BA0)($at)  
/* 002BC 80865FCC 46809120 */  cvt.s.w $f4, $f18                  
/* 002C0 80865FD0 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 002C4 80865FD4 44814000 */  mtc1    $at, $f8                   ## $f8 = 1.00
/* 002C8 80865FD8 3C014120 */  lui     $at, 0x4120                ## $at = 41200000
/* 002CC 80865FDC 44819000 */  mtc1    $at, $f18                  ## $f18 = 10.00
/* 002D0 80865FE0 3C014000 */  lui     $at, 0x4000                ## $at = 40000000
/* 002D4 80865FE4 46062002 */  mul.s   $f0, $f4, $f6              
/* 002D8 80865FE8 C4820158 */  lwc1    $f2, 0x0158($a0)           ## 00000158
/* 002DC 80865FEC 000240C0 */  sll     $t0, $v0,  3               
/* 002E0 80865FF0 01024021 */  addu    $t0, $t0, $v0              
/* 002E4 80865FF4 00084080 */  sll     $t0, $t0,  2               
/* 002E8 80865FF8 01024023 */  subu    $t0, $t0, $v0              
/* 002EC 80865FFC 2509FEE8 */  addiu   $t1, $t0, 0xFEE8           ## $t1 = FFFFFEE8
/* 002F0 80866000 46000002 */  mul.s   $f0, $f0, $f0              
/* 002F4 80866004 46004281 */  sub.s   $f10, $f8, $f0             
/* 002F8 80866008 44814000 */  mtc1    $at, $f8                   ## $f8 = 2.00
/* 002FC 8086600C 3C018086 */  lui     $at, %hi(D_80867BA4)       ## $at = 80860000
/* 00300 80866010 460C5402 */  mul.s   $f16, $f10, $f12           
/* 00304 80866014 46024281 */  sub.s   $f10, $f8, $f2             
/* 00308 80866018 46128100 */  add.s   $f4, $f16, $f18            
/* 0030C 8086601C 4600218D */  trunc.w.s $f6, $f4                   
/* 00310 80866020 44193000 */  mfc1    $t9, $f6                   
/* 00314 80866024 00000000 */  nop
/* 00318 80866028 A4990164 */  sh      $t9, 0x0164($a0)           ## 00000164
/* 0031C 8086602C C4307BA4 */  lwc1    $f16, %lo(D_80867BA4)($at) 
/* 00320 80866030 28410010 */  slti    $at, $v0, 0x0010           
/* 00324 80866034 46105482 */  mul.s   $f18, $f10, $f16           
/* 00328 80866038 46121100 */  add.s   $f4, $f2, $f18             
/* 0032C 8086603C 10200002 */  beq     $at, $zero, .L80866048     
/* 00330 80866040 E4840158 */  swc1    $f4, 0x0158($a0)           ## 00000158
/* 00334 80866044 A0890168 */  sb      $t1, 0x0168($a0)           ## 00000168
.L80866048:
/* 00338 80866048 28610009 */  slti    $at, $v1, 0x0009           
/* 0033C 8086604C 1420000E */  bne     $at, $zero, .L80866088     
/* 00340 80866050 00601025 */  or      $v0, $v1, $zero            ## $v0 = 0000FFFF
/* 00344 80866054 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 00348 80866058 44813000 */  mtc1    $at, $f6                   ## $f6 = 1.00
/* 0034C 8086605C C480015C */  lwc1    $f0, 0x015C($a0)           ## 0000015C
/* 00350 80866060 3C013E80 */  lui     $at, 0x3E80                ## $at = 3E800000
/* 00354 80866064 4606003C */  c.lt.s  $f0, $f6                   
/* 00358 80866068 00000000 */  nop
/* 0035C 8086606C 45020012 */  bc1fl   .L808660B8                 
/* 00360 80866070 28410008 */  slti    $at, $v0, 0x0008           
/* 00364 80866074 44814000 */  mtc1    $at, $f8                   ## $f8 = 0.25
/* 00368 80866078 94820166 */  lhu     $v0, 0x0166($a0)           ## 00000166
/* 0036C 8086607C 46080280 */  add.s   $f10, $f0, $f8             
/* 00370 80866080 1000000C */  beq     $zero, $zero, .L808660B4   
/* 00374 80866084 E48A015C */  swc1    $f10, 0x015C($a0)          ## 0000015C
.L80866088:
/* 00378 80866088 C480015C */  lwc1    $f0, 0x015C($a0)           ## 0000015C
/* 0037C 8086608C 44808000 */  mtc1    $zero, $f16                ## $f16 = 0.00
/* 00380 80866090 3C013E00 */  lui     $at, 0x3E00                ## $at = 3E000000
/* 00384 80866094 4600803C */  c.lt.s  $f16, $f0                  
/* 00388 80866098 00000000 */  nop
/* 0038C 8086609C 45020006 */  bc1fl   .L808660B8                 
/* 00390 808660A0 28410008 */  slti    $at, $v0, 0x0008           
/* 00394 808660A4 44819000 */  mtc1    $at, $f18                  ## $f18 = 0.12
/* 00398 808660A8 94820166 */  lhu     $v0, 0x0166($a0)           ## 00000166
/* 0039C 808660AC 46120101 */  sub.s   $f4, $f0, $f18             
/* 003A0 808660B0 E484015C */  swc1    $f4, 0x015C($a0)           ## 0000015C
.L808660B4:
/* 003A4 808660B4 28410008 */  slti    $at, $v0, 0x0008           
.L808660B8:
/* 003A8 808660B8 10200003 */  beq     $at, $zero, .L808660C8     
/* 003AC 808660BC 00000000 */  nop
/* 003B0 808660C0 A0800168 */  sb      $zero, 0x0168($a0)         ## 00000168
/* 003B4 808660C4 94820166 */  lhu     $v0, 0x0166($a0)           ## 00000166
.L808660C8:
/* 003B8 808660C8 54400004 */  bnel    $v0, $zero, .L808660DC     
/* 003BC 808660CC 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 003C0 808660D0 0C00B55C */  jal     Actor_Kill
              
/* 003C4 808660D4 A48A0166 */  sh      $t2, 0x0166($a0)           ## 00000166
/* 003C8 808660D8 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L808660DC:
/* 003CC 808660DC 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 003D0 808660E0 03E00008 */  jr      $ra                        
/* 003D4 808660E4 00000000 */  nop


