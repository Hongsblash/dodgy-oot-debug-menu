glabel ObjHana_Init
/* 00000 80B93860 27BDFFC8 */  addiu   $sp, $sp, 0xFFC8           ## $sp = FFFFFFC8
/* 00004 80B93864 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 00008 80B93868 AFB10018 */  sw      $s1, 0x0018($sp)           
/* 0000C 80B9386C AFB00014 */  sw      $s0, 0x0014($sp)           
/* 00010 80B93870 AFA5003C */  sw      $a1, 0x003C($sp)           
/* 00014 80B93874 848E001C */  lh      $t6, 0x001C($a0)           ## 0000001C
/* 00018 80B93878 3C0580B9 */  lui     $a1, %hi(D_80B93AD4)       ## $a1 = 80B90000
/* 0001C 80B9387C 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 00020 80B93880 31CF0003 */  andi    $t7, $t6, 0x0003           ## $t7 = 00000000
/* 00024 80B93884 A7AF0032 */  sh      $t7, 0x0032($sp)           
/* 00028 80B93888 0C01E037 */  jal     Actor_ProcessInitChain
              
/* 0002C 80B9388C 24A53AD4 */  addiu   $a1, $a1, %lo(D_80B93AD4)  ## $a1 = 80B93AD4
/* 00030 80B93890 87B80032 */  lh      $t8, 0x0032($sp)           
/* 00034 80B93894 3C0880B9 */  lui     $t0, %hi(D_80B93AA4)       ## $t0 = 80B90000
/* 00038 80B93898 25083AA4 */  addiu   $t0, $t0, %lo(D_80B93AA4)  ## $t0 = 80B93AA4
/* 0003C 80B9389C 0018C900 */  sll     $t9, $t8,  4               
/* 00040 80B938A0 03288821 */  addu    $s1, $t9, $t0              
/* 00044 80B938A4 8E250004 */  lw      $a1, 0x0004($s1)           ## 00000004
/* 00048 80B938A8 0C00B58B */  jal     Actor_SetScale
              
/* 0004C 80B938AC 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00050 80B938B0 C6240008 */  lwc1    $f4, 0x0008($s1)           ## 00000008
/* 00054 80B938B4 2605014C */  addiu   $a1, $s0, 0x014C           ## $a1 = 0000014C
/* 00058 80B938B8 E60400BC */  swc1    $f4, 0x00BC($s0)           ## 000000BC
/* 0005C 80B938BC 8629000C */  lh      $t1, 0x000C($s1)           ## 0000000C
/* 00060 80B938C0 8FA4003C */  lw      $a0, 0x003C($sp)           
/* 00064 80B938C4 05200015 */  bltz    $t1, .L80B9391C            
/* 00068 80B938C8 00000000 */  nop
/* 0006C 80B938CC 0C0170D9 */  jal     ActorCollider_AllocCylinder
              
/* 00070 80B938D0 AFA50024 */  sw      $a1, 0x0024($sp)           
/* 00074 80B938D4 3C0780B9 */  lui     $a3, %hi(D_80B93A70)       ## $a3 = 80B90000
/* 00078 80B938D8 24E73A70 */  addiu   $a3, $a3, %lo(D_80B93A70)  ## $a3 = 80B93A70
/* 0007C 80B938DC 8FA4003C */  lw      $a0, 0x003C($sp)           
/* 00080 80B938E0 8FA50024 */  lw      $a1, 0x0024($sp)           
/* 00084 80B938E4 0C01712B */  jal     ActorCollider_InitCylinder
              
/* 00088 80B938E8 02003025 */  or      $a2, $s0, $zero            ## $a2 = 00000000
/* 0008C 80B938EC 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00090 80B938F0 0C0189B7 */  jal     ActorCollider_Cylinder_Update
              
/* 00094 80B938F4 8FA50024 */  lw      $a1, 0x0024($sp)           
/* 00098 80B938F8 862A000C */  lh      $t2, 0x000C($s1)           ## 0000000C
/* 0009C 80B938FC 3C0680B9 */  lui     $a2, %hi(D_80B93A9C)       ## $a2 = 80B90000
/* 000A0 80B93900 24C63A9C */  addiu   $a2, $a2, %lo(D_80B93A9C)  ## $a2 = 80B93A9C
/* 000A4 80B93904 A60A018C */  sh      $t2, 0x018C($s0)           ## 0000018C
/* 000A8 80B93908 862B000E */  lh      $t3, 0x000E($s1)           ## 0000000E
/* 000AC 80B9390C 26040098 */  addiu   $a0, $s0, 0x0098           ## $a0 = 00000098
/* 000B0 80B93910 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 000B4 80B93914 0C0187B5 */  jal     func_80061ED4              
/* 000B8 80B93918 A60B018E */  sh      $t3, 0x018E($s0)           ## 0000018E
.L80B9391C:
/* 000BC 80B9391C 3C0C80B9 */  lui     $t4, %hi(D_80B93AC4)       ## $t4 = 80B90000
/* 000C0 80B93920 258C3AC4 */  addiu   $t4, $t4, %lo(D_80B93AC4)  ## $t4 = 80B93AC4
/* 000C4 80B93924 162C0007 */  bne     $s1, $t4, .L80B93944       
/* 000C8 80B93928 3C0D8016 */  lui     $t5, %hi(gSaveContext+0xedc)
/* 000CC 80B9392C 95ADF53C */  lhu     $t5, %lo(gSaveContext+0xedc)($t5)
/* 000D0 80B93930 31AE0001 */  andi    $t6, $t5, 0x0001           ## $t6 = 00000000
/* 000D4 80B93934 51C00004 */  beql    $t6, $zero, .L80B93948     
/* 000D8 80B93938 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 000DC 80B9393C 0C00B55C */  jal     Actor_Kill
              
/* 000E0 80B93940 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
.L80B93944:
/* 000E4 80B93944 8FBF001C */  lw      $ra, 0x001C($sp)           
.L80B93948:
/* 000E8 80B93948 8FB00014 */  lw      $s0, 0x0014($sp)           
/* 000EC 80B9394C 8FB10018 */  lw      $s1, 0x0018($sp)           
/* 000F0 80B93950 03E00008 */  jr      $ra                        
/* 000F4 80B93954 27BD0038 */  addiu   $sp, $sp, 0x0038           ## $sp = 00000000
