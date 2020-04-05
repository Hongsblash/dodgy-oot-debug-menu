glabel func_80904108
/* 071C8 80904108 27BDFF88 */  addiu   $sp, $sp, 0xFF88           ## $sp = FFFFFF88
/* 071CC 8090410C AFBF003C */  sw      $ra, 0x003C($sp)           
/* 071D0 80904110 AFB00038 */  sw      $s0, 0x0038($sp)           
/* 071D4 80904114 AFA40078 */  sw      $a0, 0x0078($sp)           
/* 071D8 80904118 AFA5007C */  sw      $a1, 0x007C($sp)           
/* 071DC 8090411C C4860324 */  lwc1    $f6, 0x0324($a0)           ## 00000324
/* 071E0 80904120 44802000 */  mtc1    $zero, $f4                 ## $f4 = 0.00
/* 071E4 80904124 3C068091 */  lui     $a2, %hi(D_8090D750)       ## $a2 = 80910000
/* 071E8 80904128 24C6D750 */  addiu   $a2, $a2, %lo(D_8090D750)  ## $a2 = 8090D750
/* 071EC 8090412C 4606203C */  c.lt.s  $f4, $f6                   
/* 071F0 80904130 27A40060 */  addiu   $a0, $sp, 0x0060           ## $a0 = FFFFFFE8
/* 071F4 80904134 4502007E */  bc1fl   .L80904330                 
/* 071F8 80904138 8FBF003C */  lw      $ra, 0x003C($sp)           
/* 071FC 8090413C 8CA50000 */  lw      $a1, 0x0000($a1)           ## 00000000
/* 07200 80904140 2407140B */  addiu   $a3, $zero, 0x140B         ## $a3 = 0000140B
/* 07204 80904144 0C031AB1 */  jal     Graph_OpenDisps              
/* 07208 80904148 00A08025 */  or      $s0, $a1, $zero            ## $s0 = 00000000
/* 0720C 8090414C 0C034213 */  jal     Matrix_Push              
/* 07210 80904150 00000000 */  nop
/* 07214 80904154 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 07218 80904158 3C19E700 */  lui     $t9, 0xE700                ## $t9 = E7000000
/* 0721C 8090415C 3C0CDB06 */  lui     $t4, 0xDB06                ## $t4 = DB060000
/* 07220 80904160 24580008 */  addiu   $t8, $v0, 0x0008           ## $t8 = 00000008
/* 07224 80904164 AE1802D0 */  sw      $t8, 0x02D0($s0)           ## 000002D0
/* 07228 80904168 AC590000 */  sw      $t9, 0x0000($v0)           ## 00000000
/* 0722C 8090416C AC400004 */  sw      $zero, 0x0004($v0)         ## 00000004
/* 07230 80904170 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 07234 80904174 8FAA007C */  lw      $t2, 0x007C($sp)           
/* 07238 80904178 358C0020 */  ori     $t4, $t4, 0x0020           ## $t4 = DB060020
/* 0723C 8090417C 244B0008 */  addiu   $t3, $v0, 0x0008           ## $t3 = 00000008
/* 07240 80904180 AE0B02D0 */  sw      $t3, 0x02D0($s0)           ## 000002D0
/* 07244 80904184 3C030001 */  lui     $v1, 0x0001                ## $v1 = 00010000
/* 07248 80904188 AC4C0000 */  sw      $t4, 0x0000($v0)           ## 00000000
/* 0724C 8090418C 006A1821 */  addu    $v1, $v1, $t2              
/* 07250 80904190 8C631DE4 */  lw      $v1, 0x1DE4($v1)           ## 00011DE4
/* 07254 80904194 8D440000 */  lw      $a0, 0x0000($t2)           ## 00000000
/* 07258 80904198 240C0020 */  addiu   $t4, $zero, 0x0020         ## $t4 = 00000020
/* 0725C 8090419C 00034023 */  subu    $t0, $zero, $v1            
/* 07260 809041A0 0008C040 */  sll     $t8, $t0,  1               
/* 07264 809041A4 0008C8C0 */  sll     $t9, $t0,  3               
/* 07268 809041A8 240B0020 */  addiu   $t3, $zero, 0x0020         ## $t3 = 00000020
/* 0726C 809041AC 240D0020 */  addiu   $t5, $zero, 0x0020         ## $t5 = 00000020
/* 07270 809041B0 240E0040 */  addiu   $t6, $zero, 0x0040         ## $t6 = 00000040
/* 07274 809041B4 240F0001 */  addiu   $t7, $zero, 0x0001         ## $t7 = 00000001
/* 07278 809041B8 AFAF0018 */  sw      $t7, 0x0018($sp)           
/* 0727C 809041BC AFAE0014 */  sw      $t6, 0x0014($sp)           
/* 07280 809041C0 AFAD0010 */  sw      $t5, 0x0010($sp)           
/* 07284 809041C4 AFAB0024 */  sw      $t3, 0x0024($sp)           
/* 07288 809041C8 AFB90020 */  sw      $t9, 0x0020($sp)           
/* 0728C 809041CC AFB8001C */  sw      $t8, 0x001C($sp)           
/* 07290 809041D0 AFAC0028 */  sw      $t4, 0x0028($sp)           
/* 07294 809041D4 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 07298 809041D8 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 0729C 809041DC AFA20058 */  sw      $v0, 0x0058($sp)           
/* 072A0 809041E0 0C0253D0 */  jal     Gfx_TwoTexScroll              
/* 072A4 809041E4 00603025 */  or      $a2, $v1, $zero            ## $a2 = 00010000
/* 072A8 809041E8 8FA90058 */  lw      $t1, 0x0058($sp)           
/* 072AC 809041EC 3C0EFA00 */  lui     $t6, 0xFA00                ## $t6 = FA000000
/* 072B0 809041F0 3C01FFC8 */  lui     $at, 0xFFC8                ## $at = FFC80000
/* 072B4 809041F4 AD220004 */  sw      $v0, 0x0004($t1)           ## 00000004
/* 072B8 809041F8 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 072BC 809041FC 3C19FF00 */  lui     $t9, 0xFF00                ## $t9 = FF000000
/* 072C0 80904200 37390080 */  ori     $t9, $t9, 0x0080           ## $t9 = FF000080
/* 072C4 80904204 244D0008 */  addiu   $t5, $v0, 0x0008           ## $t5 = 00000008
/* 072C8 80904208 AE0D02D0 */  sw      $t5, 0x02D0($s0)           ## 000002D0
/* 072CC 8090420C AC4E0000 */  sw      $t6, 0x0000($v0)           ## 00000000
/* 072D0 80904210 8FAF0078 */  lw      $t7, 0x0078($sp)           
/* 072D4 80904214 3C18FB00 */  lui     $t8, 0xFB00                ## $t8 = FB000000
/* 072D8 80904218 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 072DC 8090421C C5E80324 */  lwc1    $f8, 0x0324($t7)           ## 00000324
/* 072E0 80904220 4600428D */  trunc.w.s $f10, $f8                  
/* 072E4 80904224 440C5000 */  mfc1    $t4, $f10                  
/* 072E8 80904228 00000000 */  nop
/* 072EC 8090422C 318D00FF */  andi    $t5, $t4, 0x00FF           ## $t5 = 00000000
/* 072F0 80904230 01A17025 */  or      $t6, $t5, $at              ## $t6 = FFC80000
/* 072F4 80904234 AC4E0004 */  sw      $t6, 0x0004($v0)           ## 00000004
/* 072F8 80904238 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 072FC 8090423C 3C01C348 */  lui     $at, 0xC348                ## $at = C3480000
/* 07300 80904240 44816000 */  mtc1    $at, $f12                  ## $f12 = -200.00
/* 07304 80904244 244F0008 */  addiu   $t7, $v0, 0x0008           ## $t7 = 00000008
/* 07308 80904248 AE0F02D0 */  sw      $t7, 0x02D0($s0)           ## 000002D0
/* 0730C 8090424C 3C018091 */  lui     $at, %hi(D_8090DD10)       ## $at = 80910000
/* 07310 80904250 AC590004 */  sw      $t9, 0x0004($v0)           ## 00000004
/* 07314 80904254 AC580000 */  sw      $t8, 0x0000($v0)           ## 00000000
/* 07318 80904258 44066000 */  mfc1    $a2, $f12                  
/* 0731C 8090425C 0C034261 */  jal     Matrix_Translate              
/* 07320 80904260 C42EDD10 */  lwc1    $f14, %lo(D_8090DD10)($at) 
/* 07324 80904264 3C018091 */  lui     $at, %hi(D_8090DD14)       ## $at = 80910000
/* 07328 80904268 C42CDD14 */  lwc1    $f12, %lo(D_8090DD14)($at) 
/* 0732C 8090426C 3C018091 */  lui     $at, %hi(D_8090DD18)       ## $at = 80910000
/* 07330 80904270 C42EDD18 */  lwc1    $f14, %lo(D_8090DD18)($at) 
/* 07334 80904274 44066000 */  mfc1    $a2, $f12                  
/* 07338 80904278 0C0342A3 */  jal     Matrix_Scale              
/* 0733C 8090427C 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 07340 80904280 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 07344 80904284 3C0CDA38 */  lui     $t4, 0xDA38                ## $t4 = DA380000
/* 07348 80904288 358C0003 */  ori     $t4, $t4, 0x0003           ## $t4 = DA380003
/* 0734C 8090428C 244B0008 */  addiu   $t3, $v0, 0x0008           ## $t3 = 00000008
/* 07350 80904290 AE0B02D0 */  sw      $t3, 0x02D0($s0)           ## 000002D0
/* 07354 80904294 AC4C0000 */  sw      $t4, 0x0000($v0)           ## 00000000
/* 07358 80904298 8FAD007C */  lw      $t5, 0x007C($sp)           
/* 0735C 8090429C 3C058091 */  lui     $a1, %hi(D_8090D764)       ## $a1 = 80910000
/* 07360 809042A0 24A5D764 */  addiu   $a1, $a1, %lo(D_8090D764)  ## $a1 = 8090D764
/* 07364 809042A4 8DA40000 */  lw      $a0, 0x0000($t5)           ## 00000000
/* 07368 809042A8 2406143F */  addiu   $a2, $zero, 0x143F         ## $a2 = 0000143F
/* 0736C 809042AC 0C0346A2 */  jal     Matrix_NewMtx              
/* 07370 809042B0 AFA2004C */  sw      $v0, 0x004C($sp)           
/* 07374 809042B4 8FA3004C */  lw      $v1, 0x004C($sp)           
/* 07378 809042B8 3C048091 */  lui     $a0, %hi(D_8090B100)       ## $a0 = 80910000
/* 0737C 809042BC 2484B100 */  addiu   $a0, $a0, %lo(D_8090B100)  ## $a0 = 8090B100
/* 07380 809042C0 AC620004 */  sw      $v0, 0x0004($v1)           ## 00000004
/* 07384 809042C4 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 07388 809042C8 0004C100 */  sll     $t8, $a0,  4               
/* 0738C 809042CC 0018CF02 */  srl     $t9, $t8, 28               
/* 07390 809042D0 244E0008 */  addiu   $t6, $v0, 0x0008           ## $t6 = 00000008
/* 07394 809042D4 AE0E02D0 */  sw      $t6, 0x02D0($s0)           ## 000002D0
/* 07398 809042D8 00195880 */  sll     $t3, $t9,  2               
/* 0739C 809042DC 3C0FDE00 */  lui     $t7, 0xDE00                ## $t7 = DE000000
/* 073A0 809042E0 3C0C8016 */  lui     $t4, 0x8016                ## $t4 = 80160000
/* 073A4 809042E4 018B6021 */  addu    $t4, $t4, $t3              
/* 073A8 809042E8 3C0100FF */  lui     $at, 0x00FF                ## $at = 00FF0000
/* 073AC 809042EC AC4F0000 */  sw      $t7, 0x0000($v0)           ## 00000000
/* 073B0 809042F0 8D8C6FA8 */  lw      $t4, 0x6FA8($t4)           ## 80166FA8
/* 073B4 809042F4 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = 00FFFFFF
/* 073B8 809042F8 00816824 */  and     $t5, $a0, $at              
/* 073BC 809042FC 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 073C0 80904300 018D7021 */  addu    $t6, $t4, $t5              
/* 073C4 80904304 01C17821 */  addu    $t7, $t6, $at              
/* 073C8 80904308 0C034221 */  jal     Matrix_Pull              
/* 073CC 8090430C AC4F0004 */  sw      $t7, 0x0004($v0)           ## 00000004
/* 073D0 80904310 8FB8007C */  lw      $t8, 0x007C($sp)           
/* 073D4 80904314 3C068091 */  lui     $a2, %hi(D_8090D778)       ## $a2 = 80910000
/* 073D8 80904318 24C6D778 */  addiu   $a2, $a2, %lo(D_8090D778)  ## $a2 = 8090D778
/* 073DC 8090431C 27A40060 */  addiu   $a0, $sp, 0x0060           ## $a0 = FFFFFFE8
/* 073E0 80904320 24071442 */  addiu   $a3, $zero, 0x1442         ## $a3 = 00001442
/* 073E4 80904324 0C031AD5 */  jal     Graph_CloseDisps              
/* 073E8 80904328 8F050000 */  lw      $a1, 0x0000($t8)           ## 00000000
/* 073EC 8090432C 8FBF003C */  lw      $ra, 0x003C($sp)           
.L80904330:
/* 073F0 80904330 8FB00038 */  lw      $s0, 0x0038($sp)           
/* 073F4 80904334 27BD0078 */  addiu   $sp, $sp, 0x0078           ## $sp = 00000000
/* 073F8 80904338 03E00008 */  jr      $ra                        
/* 073FC 8090433C 00000000 */  nop


