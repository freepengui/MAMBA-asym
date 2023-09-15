/*
Copyright (c) 2017, The University of Bristol, Senate House, Tyndall Avenue, Bristol, BS8 1TH, United Kingdom.
Copyright (c) 2021, COSIC-KU Leuven, Kasteelpark Arenberg 10, bus 2452, B-3001 Leuven-Heverlee, Belgium.

All rights reserved
*/

#include "MASCOTTriples.h"
#include "LSSS/PRSS.h"
#include "Math/gfp.h"
#include "OT/SimpleROT.h"
#include "Tools/Crypto.h"
#include "config.h"
#include <unistd.h>

using namespace std::chrono;

MASCOTTriples::MASCOTTriples(Player &Pl, bigint fsize, gfp delta, unsigned int connection) : P(Pl)
{
  last_thread= 0;
  myThreadNb= connection;
  totalTriplesGenerated= 0;
  size_field= fsize;
  bit_size_field= numBits(size_field);
  keyGen_sec= 128;
  gf2n_to_gfp_sec= ceil((bit_size_field + keyGen_sec) / (float) OT_comp_sec);
  tau= 1;
  cout << "[MASCOTTriples - constructor] tau (number of triples to be generated per output triples) is set to " << tau << endl;
  nbTriples= 1 << 9;
  cout << "[MASCOTTriples - constructor] Batch size for triple generation set to " << nbTriples << endl;
  gf2n::init_field(OT_comp_sec);
  gfp::init_field(size_field);

  mac_usable_triples.resize(3);
  usable_triples.resize(3);

  nbPlayers= P.nplayers();
  nbOTs= (unsigned long int) nbTriples * (unsigned long int) gf2n_to_gfp_sec * (unsigned long int) tau * (unsigned long int) bit_size_field;

  nbOTsOutput= tau * bit_size_field * nbTriples;
  // Initalize every vector once and for all (trying to solve memory issue)
  triples.resize(3);
  triples[0].resize(tau * nbTriples);
  triples[1].resize(nbTriples);
  triples[2].resize(tau * nbTriples);

  ROTS.resize(nbPlayers);
  ROTR.resize(nbPlayers);
  COPE_S.resize(nbPlayers);
  COPE_R.resize(nbPlayers);

  g.resize(bit_size_field);
  for (unsigned int i= 0; i < bit_size_field; i++)
    {
      g[i].assign(2);
      g[i].power(i);
    }

  combined_triples.resize(5);
  mac_combined_triples.resize(5);
  for (unsigned int i= 0; i < 5; i++)
    {
      if (i != 1)
        {
          combined_triples[i].resize(nbTriples);
        }
      mac_combined_triples[i].resize(nbTriples);
    }

  // Initialize every PRG
  G.ReSeed(myThreadNb + 1);

  // G.ReSeed(0); // PRG seed

  uint8_t seed[SEED_SIZE];
  AgreeRandom(P, seed, SEED_SIZE, 1);
  G2.SetSeedFromRandom(seed);

  isInit= false;
  Delta= delta;
}

MASCOTTriples::MASCOTTriples(Player &Pl, bigint fsize, unsigned int connection) : P(Pl)
{
  last_thread= 0;
  myThreadNb= connection;
  totalTriplesGenerated= 0;
  size_field= fsize;                   // 域size 大小
  bit_size_field= numBits(size_field); // 比特数 120
  keyGen_sec= 128;
  gf2n_to_gfp_sec= ceil((bit_size_field + keyGen_sec) / (float) OT_comp_sec); // 这个表示？
  // tau= 3;
  tau= 1; // 非对称
  cout << "[MASCOTTriples - constructor] tau (number of triples to be generated per output triples) is set to " << tau << endl;
  nbTriples= 1 << 9; // 64
  cout << "[MASCOTTriples - constructor] Batch size for triple generation set to " << nbTriples << endl;
  gf2n::init_field(OT_comp_sec);
  gfp::init_field(size_field);

  mac_usable_triples.resize(3);
  usable_triples.resize(3);
  Delta.assign(0); // 这里不一样

  nbPlayers= P.nplayers(); // 2
  nbOTs= (unsigned long int) nbTriples * (unsigned long int) gf2n_to_gfp_sec * (unsigned long int) tau * (unsigned long int) bit_size_field;
  // 需要的OTs 数量吧？ 46080    = 64 * 2 * 3 * 120

  nbOTsOutput= tau * bit_size_field * nbTriples; // OT 输出数量， 三元组数量*tau*域的比特数 = nbOTs / 2 = 23040
  // Initalize every vector once and for all (trying to solve memory issue)
  triples.resize(3);                  // a  b  c
  triples[0].resize(tau * nbTriples); // a[tau-tau-...-tau]
  triples[1].resize(nbTriples);       // b[1-1-..-1]
  triples[2].resize(tau * nbTriples); // c[tau-tau-...-tau]

  ROTS.resize(nbPlayers); // ROT sender , 2,两个方向
  ROTR.resize(nbPlayers);
  COPE_S.resize(nbPlayers);
  COPE_R.resize(nbPlayers);

  g.resize(bit_size_field); // vector<gfp> g; 比特分解相关？
  for (unsigned int i= 0; i < bit_size_field; i++)
    {
      g[i].assign(2); // base 2
      g[i].power(i);  // 幂
    }

  combined_triples.resize(5);     // 组合后 5 个 a c b a^ c^
  mac_combined_triples.resize(5); // 带上mac
  for (unsigned int i= 0; i < 5; i++)
    {
      if (i != 1)
        {
          combined_triples[i].resize(nbTriples);
        }
      mac_combined_triples[i].resize(nbTriples);
    }

  // Initialize every PRG
  G.ReSeed(0); // PRG seed

  uint8_t seed[SEED_SIZE];            // SEED_SIZE AES_BLK_SIZE 16
  AgreeRandom(P, seed, SEED_SIZE, 1); // 协商随机数?
  G2.SetSeedFromRandom(seed);         // 从随机数生成种子

  isInit= false; // 是否初始化
}

// Player 0 says if we should stop triple generation
int check_exit(const Player &P, int myThreadNb, bool kg)
{
  int result= 0;
  string ss= "-";
  if (P.whoami() == 0)
    {
      if (not kg)
        {
          result= 1;
          ss= "E";
        }
      P.send_to_player(1,ss,1);
    }
  else
    {
      P.receive_from_player(0, ss, 1);
      if (ss.compare("E") == 0)
        {
          result= 1;
        }
    }
  return result;
}

void MASCOTTriples::execute(std::mutex &mtx, vector<list<gfp>> &u_triples, vector<list<gfp>> &mac_u_triples, bool &keepGoing, unsigned long int *totalProduced)
{
  gfp::init_field(size_field);
  cout << "Starting producing in THREAD" << myThreadNb << endl;
  // bool kg = true;
  unsigned int localProduced= 0;
  while (true)
    {
      if (check_exit(P, myThreadNb, keepGoing) == 1)
        {
          break;
        }
      Multiply();
      Combine();
      Authenticate();
      Sacrifice();
      mtx.lock();
      *totalProduced= *totalProduced + nbTriples;
      localProduced+= nbTriples;
      for (unsigned int j= 0; j < 3; j++)
        {
          while (!usable_triples[j].empty())
            {
              u_triples[j].push_back(usable_triples[j].back());
              usable_triples[j].pop_back();

              mac_u_triples[j].push_back(mac_usable_triples[j].back());
              mac_usable_triples[j].pop_back();
            }
        }
      if (*totalProduced % 100000 < nbTriples)
        {
          cout << "Produced " << *totalProduced << " triples so far " << endl;
        }
      mtx.unlock();
    }

  cout << "Ending production in THREAD" << myThreadNb << endl;
  cout << "Produced " << *totalProduced << " triples in total in this field" << endl;
}

void MASCOTTriples::Multiply()
{
  // 1st step of MASCOT
  for (unsigned int i= 0; i < nbTriples; i++)
    {
      triples[1][i].randomize(G); // 三元组b随机生成
      for (unsigned int j= 0; j < tau; j++)
        {
          triples[0][i * tau + j].randomize(G); // 三元组 a 随机生成 向量（共tau*ntriples, 每 tau=3 个属于一份三元组）
        }
    }

  // 2nd step of MASCOT
  unsigned int whoami= P.whoami();
  CryptoPP::RandomPool RNG;

  // Extract nbOTs bits from triples[0] (which represents a)   // a比特分解 ？
  //  bigint tmp;
  //  unsigned long int currIndex;    // 代表什么意思呢？
  //  BitVector bit_triples0(nbOTs);  // 比特分解后 比特元组 nbOTs = nbTriples * tau * bit_size_field * gf2n_to_gfp_sec
  //  nbOTs个bits,
  //  class BitVector  {uint8_t *bytes; size_t nbytes;  size_t nbits; size_t length;} 每个OT对应一个BitVector
  //  for (unsigned long int i= 0; i < nbTriples * tau; i++)
  //    {
  //      vector<int> tmp_bit_array(bit_size_field);          // 临时的bit数组，大小为 bit_size_field=120 这个120 到底是什么意思，为什么要120，从最初的fsize1来的
  //      to_bigint(tmp, triples[0][i], false);               // gfp -> bigint
  //      bigint_to_bits(tmp_bit_array, tmp);                 // bigint -> bits  tem_bit_array中是120个比特（int 0/1）
  //      for (unsigned int j= 0; j < bit_size_field; j++)
  //        {
  //          for (unsigned int jj= 0; jj < gf2n_to_gfp_sec; jj++)  // g2fn -> gfp  = 2
  //            {
  //              currIndex= (i * bit_size_field * gf2n_to_gfp_sec) + (j * gf2n_to_gfp_sec) + jj;
  //              bit_triples0.set_bit(currIndex, tmp_bit_array[j]);              // 比特分解后 存储
  //            }
  //        }
  //    }

  // class BitVector  {uint8_t *bytes; size_t nbytes;  size_t nbits; size_t length;} 每个OT对应一个BitVector
  bigint tmp;
  unsigned long int currIndex;     // 代表什么意思呢？
  BitVector bit_triples_b1(nbOTs); // 对b1比特分解后 比特元组 nbOTs = nbTriples * tau * bit_size_field * gf2n_to_gfp_sec | nbOTs个bits,
  BitVector bit_triples_a1(nbOTs); // 对a1比特分解
  // 非对称 a0*b1 + b0*a1 ,所以 b1 a1都要比特分解，两次ROT中作为选择比特，只有序号为1的参与方需要分解，所以需要判断

  if (whoami == 1)
    {
      // 对b1比特分解
      for (unsigned long int i= 0; i < nbTriples; i++)
        {
          vector<int> tmp_bit_array(bit_size_field); // 临时的bit数组，大小为 bit_size_field=120 这个120 到底是什么意思，为什么要120，从最初的fsize1来的
          to_bigint(tmp, triples[1][i], false);      // gfp -> bigint
          bigint_to_bits(tmp_bit_array, tmp);        // bigint -> bits  tem_bit_array中是120个比特（int 0/1）
          for (unsigned int j= 0; j < bit_size_field; j++)
            {
              for (unsigned int jj= 0; jj < gf2n_to_gfp_sec; jj++) // g2fn -> gfp  = 2
                {
                  currIndex= (i * bit_size_field * gf2n_to_gfp_sec) + (j * gf2n_to_gfp_sec) + jj;
                  bit_triples_b1.set_bit(currIndex, tmp_bit_array[j]); // 比特分解后 存储
                }
            }
        }
      // a1比特分解
      for (unsigned long int i= 0; i < nbTriples * tau; i++)
        {
          vector<int> tmp_bit_array(bit_size_field); // 临时的bit数组，大小为 bit_size_field=120 这个120 到底是什么意思，为什么要120，从最初的fsize1来的
          to_bigint(tmp, triples[0][i], false);      // gfp -> bigint
          bigint_to_bits(tmp_bit_array, tmp);        // bigint -> bits  tem_bit_array中是120个比特（int 0/1）
          for (unsigned int j= 0; j < bit_size_field; j++)
            {
              for (unsigned int jj= 0; jj < gf2n_to_gfp_sec; jj++) // gf2n_to_gfp_sec = 2
                {
                  currIndex= (i * bit_size_field * gf2n_to_gfp_sec) + (j * gf2n_to_gfp_sec) + jj;
                  bit_triples_a1.set_bit(currIndex, tmp_bit_array[j]); // 比特分解后 存储
                }
            }
        }
    }

  /* Create choicebits for Delta */ // Delta 是什么？ 用于干什么BaseOT
  if (!isInit)
    {
      vector<int> choicebits(OT_comp_sec);          // 选择比特vector 跟OT的计算安全参数有关？,128 ,
      for (unsigned int i= 0; i < OT_comp_sec; i++) //
        {
          choicebits[i]= G.get_uchar() & 1; // 每个选择比特 设置随机的无符号char 保留最后一位
        }
      /*(a) and (b), ROT stuff */
      // We do first the init COT between i and whoami
      //  for (unsigned int i= 0; i < nbPlayers; i++)
      //    {
      //      if (i < whoami)                             // 序号我
      //        {
      //          cout << "1 init" << endl;
      //          ROTS[i].init(P, i, RNG, choicebits, 1); // 我作为ROT的Sender 初始化 receiver是 i   CryptoPP::RandomPool RNG 随机数
      //        }
      //      else if (i > whoami)                        // 序号大于本参与方
      //        {
      //          cout << "2 init" << endl;
      //          ROTR[i].init(P, i, RNG, 1);             // ROT 作为 Receiver 初始化
      //        }
      //    }

      // // The other way around                          // 相反反向的OT
      // for (unsigned int i= 0; i < nbPlayers; i++)
      //   {
      //     if (i > whoami)
      //       {
      //         cout << "3 init" << endl;
      //         ROTS[i].init(P, i, RNG, choicebits, 1); // 大于本参与 作为Sender
      //       }
      //     else if (i < whoami)
      //       {
      //         cout << "4 init" << endl;
      //         ROTR[i].init(P, i, RNG, 1);             // 小于参与方 作为Receiver
      //       }
      //   }

      // 非对称
      //  如果是0(半诚实方)，作为ROT的发送者
      if (whoami == 0)
        {
          // cout << "init Base OT xxx" << endl;
          ROTS[0].init(P, 1, RNG, choicebits, 1); // 半诚实0方作为Sender 初始化 CryptoPP::RandomPool RNG 随机数
          // cout << " init 一个好了 " << endl;
          ROTS[1].init(P, 1, RNG, choicebits, 1); // 半诚实0方作为Sender 初始化
          // cout << "init close" << endl;
        }
      else
        {
          cout << "Base OT" << endl;
          ROTR[0].init(P, 0, RNG, 1); // 恶意方1 作为ROT Receiver 初始化
          cout << "ROTR Base finish one " << endl;
          ROTR[1].init(P, 0, RNG, 1); // 恶意方1 作为ROT Receiver 初始化
          cout << "Base OT over" << endl;
        }
    }

  /*Doing calls to OT next round and converting output from gf2n to gfp*/
  vector<vector<gfp>> dh_receiver(nbPlayers, vector<gfp>(nbOTsOutput));                                        // dh_receiver? 这不是跟 ROT 的sender是同一个参与方吗？
  vector<vector<vector<gfp>>> out_vec_sender_gfp(nbPlayers, vector<vector<gfp>>(nbOTsOutput, vector<gfp>(2))); // 输出向量 sender
  {
    vector<vector<gf2n>> out_vec_receiver(nbPlayers, vector<gf2n>(nbOTs)); // ROT receiver 的输出 gf2n
    // 两个参与方，每个参与方nbOTs, 每次OT sender两个gf2n的输出
    vector<vector<vector<gf2n>>> out_vec_sender(nbPlayers, vector<vector<gf2n>>(nbOTs, vector<gf2n>(2)));

    vector<gf2n> out_vec_receiver_D(2 * nbOTs); // ROT receiver 的输出 gf2n
    vector<vector<gf2n>> out_vec_sender_D(2 * nbOTs, vector<gf2n>(2));

    // for (unsigned int i= 0; i < nbPlayers; i++)
    //   {
    //     if (i < whoami)
    //       {
    //         cout << i << "iteration1" << endl;
    //         ROTR[i].next_iteration(P, nbOTs, bit_triples0, out_vec_receiver[i]); // 这里把a 比特分解后的向量传入了，receiver输出
    //         cout << "log" << endl;
    //       }
    //     else if (i > whoami)
    //       {
    //         cout << i << "iteration2" << endl;
    //         ROTS[i].next_iteration(P, nbOTs, out_vec_sender[i]);  // ROT 的sender 的输出
    //       }
    //   }

    // for (unsigned int i= 0; i < nbPlayers; i++)
    //   {
    //     if (i > whoami)
    //       {
    //         cout << i << "iteration3" << endl;
    //         ROTR[i].next_iteration(P, nbOTs, bit_triples0, out_vec_receiver[i]);
    //       }
    //     else if (i < whoami)
    //       {
    //         cout << i << "iteration4" << endl;
    //         ROTS[i].next_iteration(P, nbOTs, out_vec_sender[i]);
    //       }
    //   }
    // 得到gf2n的输出，下面将 gf2n 转 gfp;  ROT 输出？

    // 非对称两方
    if (whoami == 0)
      {
        // cout << "执行ROT Send" << endl;
        ROTS[0].next_iteration(P, nbOTs, out_vec_sender[0]); // ROT 的sender 的输出
        ROTS[1].next_iteration(P, nbOTs, out_vec_sender[1]); // ROT 的sender 的输出   //Sender: q_h_0  q_h_1
        // ROTS[1].next_iteration(P, 2 * nbOTs, out_vec_sender_D); // ROT 的sender 的输出   //Sender: q_h_0  q_h_1
        // for (unsigned int i= 0; i < nbOTs; i++)
        //   {
        //     out_vec_sender[0][i]= out_vec_sender_D[i];
        //     out_vec_sender[1][i]= out_vec_sender_D[i + nbOTs];
        //   }

        // cout << "执行完 Send" << endl;
      }
    else
      {
        // cout << "执行 ROT Receive" << endl;
        

        ROTR[0].next_iteration(P, nbOTs, bit_triples_b1, out_vec_receiver[0]); // 这里把b1 比特分解后的向量传入了，receiver输出
        ROTR[1].next_iteration(P, nbOTs, bit_triples_a1, out_vec_receiver[1]); // Receiver :  q_h_b(h)
        // ROTR[0].next_iteration(P, 2 * nbOTs, bit_triples, out_vec_receiver_D); // 这里把b1 比特分解后的向量传入了，receiver输出

        // cout << "执行完 Receive" << endl;
      }

    // Converting gf2n output of OTs to gfp        // OT的输出 转换 ，GF(2^n) - > GF(p)
    uint8_t *buffer_rcv;  // 字符数组（字符串）
    uint8_t *buffer_rcv1; // 字符数组（字符串）

    // uint8_t *buffer_send0;
    // uint8_t *buffer_send1;
    uint8_t *buffer_send00;
    uint8_t *buffer_send01;
    uint8_t *buffer_send10;
    uint8_t *buffer_send11;
    buffer_rcv= (uint8_t *) malloc(16 * gf2n_to_gfp_sec * sizeof(uint8_t)); // 申请空间，gf2n_to_gfp_sec = 2  16是怎么来的？
    buffer_send00= (uint8_t *) malloc(16 * gf2n_to_gfp_sec * sizeof(uint8_t));
    buffer_send01= (uint8_t *) malloc(16 * gf2n_to_gfp_sec * sizeof(uint8_t));

    buffer_rcv1= (uint8_t *) malloc(16 * gf2n_to_gfp_sec * sizeof(uint8_t)); // 申请空间，gf2n_to_gfp_sec = 2  16是怎么来的？
    buffer_send10= (uint8_t *) malloc(16 * gf2n_to_gfp_sec * sizeof(uint8_t));
    buffer_send11= (uint8_t *) malloc(16 * gf2n_to_gfp_sec * sizeof(uint8_t));

    // for (unsigned int i= 0; i < nbPlayers; i++)
    //   {
    //     if (i != P.whoami())
    //       {
    //         for (unsigned int j= 0; j < nbOTsOutput; j++)
    //           {
    //             for (unsigned int jj= 0; jj < gf2n_to_gfp_sec; jj++)
    //               {
    //                 currIndex= j * gf2n_to_gfp_sec + jj;
    //                 out_vec_receiver[i][currIndex].store_into_buffer(&buffer_rcv[jj * 16]);
    //                 out_vec_sender[i][currIndex][0].store_into_buffer(&buffer_send0[jj * 16]);
    //                 out_vec_sender[i][currIndex][1].store_into_buffer(&buffer_send1[jj * 16]);
    //               }
    //             bigintFromBytes(tmp, buffer_rcv, 16 * gf2n_to_gfp_sec);
    //             dh_receiver[i][j].assign(tmp);

    //             bigintFromBytes(tmp, buffer_send0, 16 * gf2n_to_gfp_sec);
    //             out_vec_sender_gfp[i][j][0].assign(tmp);

    //             bigintFromBytes(tmp, buffer_send1, 16 * gf2n_to_gfp_sec);
    //             out_vec_sender_gfp[i][j][1].assign(tmp);
    //           }
    //       }
    //   }

    if (whoami == 0)
      {
        for (unsigned int j= 0; j < nbOTsOutput; j++)
          {
            for (unsigned int jj= 0; jj < gf2n_to_gfp_sec; jj++)
              {
                currIndex= j * gf2n_to_gfp_sec + jj;

                out_vec_sender[0][currIndex][0].store_into_buffer(&buffer_send00[jj * 16]); // q_h_0
                out_vec_sender[0][currIndex][1].store_into_buffer(&buffer_send01[jj * 16]); // q_h_1

                // 第2批OT
                out_vec_sender[1][currIndex][0].store_into_buffer(&buffer_send10[jj * 16]);
                out_vec_sender[1][currIndex][1].store_into_buffer(&buffer_send11[jj * 16]);
              }
            bigintFromBytes(tmp, buffer_send00, 16 * gf2n_to_gfp_sec);
            out_vec_sender_gfp[0][j][0].assign(tmp);
            bigintFromBytes(tmp, buffer_send01, 16 * gf2n_to_gfp_sec);
            out_vec_sender_gfp[0][j][1].assign(tmp);

            bigintFromBytes(tmp, buffer_send10, 16 * gf2n_to_gfp_sec);
            out_vec_sender_gfp[1][j][0].assign(tmp);
            bigintFromBytes(tmp, buffer_send11, 16 * gf2n_to_gfp_sec);
            out_vec_sender_gfp[1][j][1].assign(tmp);
          }
      }
    else
      {
        for (unsigned int j= 0; j < nbOTsOutput; j++)
          {
            for (unsigned int jj= 0; jj < gf2n_to_gfp_sec; jj++)
              {
                currIndex= j * gf2n_to_gfp_sec + jj;
                out_vec_receiver[0][currIndex].store_into_buffer(&buffer_rcv[jj * 16]);
                out_vec_receiver[1][currIndex].store_into_buffer(&buffer_rcv1[jj * 16]);
              }
            bigintFromBytes(tmp, buffer_rcv, 16 * gf2n_to_gfp_sec);
            dh_receiver[0][j].assign(tmp); // dh_receiver = q_h_b(h)

            bigintFromBytes(tmp, buffer_rcv1, 16 * gf2n_to_gfp_sec);
            dh_receiver[1][j].assign(tmp);
          }
      }

    free(buffer_rcv);
    free(buffer_send00);
    free(buffer_send01);

    free(buffer_rcv1);
    free(buffer_send10);
    free(buffer_send11);
  }

  /* (c) - compute the dh and send them */ // dh = q_0,h - q_1,h + bj
  unsigned int indexRound;                 //
  int currBit;                             //
  gfp tmp_gfp1, tmp_gfp2;
  // 一个方向上
  // for (unsigned int i= 0; i < nbPlayers; i++)
  //   {
  //     if (i < P.whoami()) // 小于本参与方
  //       {
  //         ostringstream os;
  //         for (unsigned int j= 0; j < nbOTsOutput; j++) // 184320
  //           {
  //             indexRound= floor(j / (tau * bit_size_field)); // j/(3*120) = 512
  //             tmp_gfp= (out_vec_sender_gfp[i][j][0] - out_vec_sender_gfp[i][j][1]) + triples[1][indexRound];
  //             //  dh =                   q_h,0      -                     q_h,1    + b
  //             tmp_gfp.output(os, false);
  //           }
  //         P.send_to_player(i, os.str(), 1);
  //       }
  //     else if (i > P.whoami()) // 大于本参与方序号
  //       {
  //         string ss;
  //         P.receive_from_player(i, ss, 1);
  //         istringstream is(ss);

  //         for (unsigned int j= 0; j < nbOTsOutput; j++)
  //           {
  //             tmp_gfp.input(is, false);
  //             currBit= bit_triples0.get_bit(j * gf2n_to_gfp_sec);
  //             if (currBit != 0)
  //               {
  //                 dh_receiver[i][j]+= tmp_gfp;
  //               }
  //           }
  //       }
  //   }

  // for (unsigned int i= 0; i < nbPlayers; i++) // 另一个方向上
  //   {
  //     if (i > P.whoami()) // dh 发送方
  //       {

  //         ostringstream os;
  //         for (unsigned int j= 0; j < nbOTsOutput; j++)
  //           {
  //             indexRound= floor(j / (tau * bit_size_field));
  //             tmp_gfp= (out_vec_sender_gfp[i][j][0] - out_vec_sender_gfp[i][j][1]) + triples[1][indexRound];
  //             tmp_gfp.output(os, false);
  //           }
  //         P.send_to_player(i, os.str(), 1);
  //       }
  //     else if (i < P.whoami()) // dh接收方
  //       {
  //         string ss;
  //         P.receive_from_player(i, ss, 1);
  //         istringstream is(ss);

  //         for (unsigned int j= 0; j < nbOTsOutput; j++)
  //           {
  //             tmp_gfp.input(is, false);
  //             currBit= bit_triples0.get_bit(j * gf2n_to_gfp_sec);
  //             if (currBit != 0)
  //               {
  //                 dh_receiver[i][j]+= tmp_gfp;
  //               }
  //           }
  //       }
  //   }

  // 非对称 计算 dh

  if (P.whoami() == 0)
    {
      ostringstream os1, os2;
      for (unsigned int j= 0; j < nbOTsOutput; j++) //
        {
          indexRound= floor(j / (tau * bit_size_field)); // j/(3*120) = 512
          tmp_gfp1= (out_vec_sender_gfp[0][j][0] - out_vec_sender_gfp[0][j][1]) + triples[0][indexRound];
          //  dh =                   q_h,0      -                     q_h,1    + a_0
          // if (j < 5)
          //   {
          //     cout << "tmp_gfp1_dh_send = " << tmp_gfp1 << endl;
          //   }

          tmp_gfp1.output(os1, false);

          tmp_gfp2= (out_vec_sender_gfp[1][j][0] - out_vec_sender_gfp[1][j][1]) + triples[1][indexRound];
          //  dh =                   q_h,0      -                     q_h,1    + b_0
          // if (j < 5)
          //   {
          //     cout << "tmp_gfp2_dh_send = " << tmp_gfp2 << endl;
          //   }

          tmp_gfp2.output(os2, false);
        }
      P.send_to_player(1, os1.str(), 1);
      P.send_to_player(1, os2.str(), 1);
    }
  else
    {
      string ss1, ss2;
      P.receive_from_player(0, ss1, 1);
      P.receive_from_player(0, ss2, 1);

      istringstream is1(ss1);
      istringstream is2(ss2);
      for (unsigned int j= 0; j < nbOTsOutput; j++)
        {
          // if (j < 5)
          //   {
          //     cout << "tmp_gfp1_dh_rec1 = " << tmp_gfp1 << endl;
          //   }

          tmp_gfp1.input(is1, false);
          // if (j < 5)
          //   {
          //     cout << "tmp_gfp1_dh_rec2 = " << tmp_gfp1 << endl;
          //   }

          currBit= bit_triples_b1.get_bit(j * gf2n_to_gfp_sec); //  currBit = b(h)
          if (currBit != 0)
            {
              dh_receiver[0][j]+= tmp_gfp1; // tmp_gfp1=dh;  t_h = q_h_b(h) + b(h) * dh  此时，dh_receive中存方 t_h
            }

          // if (j < 5)
          //   {
          //     cout << "tmp_gfp1_dh_rec3 = " << tmp_gfp2 << endl;
          //   }
          tmp_gfp2.input(is2, false);
          // if (j < 5)
          //   {
          //     cout << "tmp_gfp1_dh_rec4 = " << tmp_gfp2 << endl;
          //   }

          currBit= bit_triples_a1.get_bit(j * gf2n_to_gfp_sec);
          if (currBit != 0)
            {
              dh_receiver[1][j]+= tmp_gfp2; // dh_receive = t_h
            }
        }
    }

  /*(3) compute c from intermediate c vectors */
  for (unsigned int i= 0; i < nbTriples * tau; i++)
    {
      currIndex= floor(i / tau);
      triples[2][i]= triples[0][i] * triples[1][currIndex];
      //         c =            a  * b  本地的
      // for (unsigned int j= 0; j < nbPlayers; j++)
      //   {
      // if (j != P.whoami())
      //   {
      //     for (unsigned int k= 0; k < bit_size_field; k++)
      //       {
      //         triples[2][i]+= g[k] *
      //                         (dh_receiver[j][k + ((i % tau) + currIndex * tau) * bit_size_field] -
      //                          out_vec_sender_gfp[j][k + ((i % tau) + currIndex * tau) * bit_size_field][0]); // 这里也要改

      //         // c = a * b(本地) + gk *(t_h - q_h_0)
      //       }
      //   }

      // 非对称
      for (unsigned int k= 0; k < bit_size_field; k++)
        {
          if (P.whoami() == 0)
            {
              triples[2][i]-= (g[k] *
                               (out_vec_sender_gfp[0][k + ((i % tau) + currIndex * tau) * bit_size_field][0] +
                                out_vec_sender_gfp[1][k + ((i % tau) + currIndex * tau) * bit_size_field][0])); // 这里也要改

              // c = a * b(本地) + gk *( - q_h_0) +gk * (- q_h_0)
            }
          else
            {
              triples[2][i]+= g[k] *
                              (dh_receiver[0][k + ((i % tau) + currIndex * tau) * bit_size_field] +
                               dh_receiver[1][k + ((i % tau) + currIndex * tau) * bit_size_field]); // 这里也要改

              // c = a * b(本地) + gk *(t_h + t_h)
            }
        }
      // cout << "i = " << i << " triples_c = " << triples[2][i] << endl;
      // }
    }
  // vector<vector<gfp>> triples_;
  // triples_.resize(3);                  // a  b  c
  // triples_[0].resize(tau * nbTriples); // a[tau-tau-...-tau]
  // triples_[1].resize(nbTriples);       // b[1-1-..-1]
  // triples_[2].resize(tau * nbTriples); // c[tau-tau-...-tau]
  // gfp tmp_gfp_a, tmp_gfp_b, tmp_gfp_c;
  // gfp a_b, c;
  // ostringstream osa, osb, osc;
  // if (P.whoami() == 0)
  //   {
  //     for (unsigned int i= 0; i < nbTriples; i++)
  //       {
  //         // indexRound= floor(j / (tau * bit_size_field)); // j/(3*120) = 512
  //         tmp_gfp_a= triples[0][i];
  //         tmp_gfp_a.output(osa, false);
  //         cout << " P0 triples_a = " << tmp_gfp_a << endl;

  //         tmp_gfp_b= triples[1][i];
  //         tmp_gfp_b.output(osb, false);
  //         cout << " P0 triples_b = " << tmp_gfp_b << endl;

  //         tmp_gfp_c= triples[2][i];
  //         tmp_gfp_c.output(osc, false);
  //         cout << " P0 triples_c = " << tmp_gfp_c << endl;
  //       }
  //     P.send_to_player(1, osa.str(), 1);
  //     // cout << "????  1" << endl;
  //     P.send_to_player(1, osb.str(), 2);
  //     // cout << "????  2" << endl;
  //     P.send_to_player(1, osc.str(), 1);
  //     // cout << "????  3" << endl;
  //   }
  // else
  //   {
  //     string ssa, ssb, ssc;
  //     P.receive_from_player(0, ssa, 1);
  //     P.receive_from_player(0, ssb, 2);
  //     P.receive_from_player(0, ssc, 1);

  //     istringstream isa(ssa);
  //     istringstream isb(ssb);
  //     istringstream isc(ssc);
  //     for (unsigned int j= 0; j < nbTriples; j++)
  //       {
  //         tmp_gfp_a.input(isa, false);
  //         triples_[0][j]= tmp_gfp_a;
  //         // cout << "P0 triples_a_ = " << tmp_gfp_a << endl;
  //         // cout << "P1 triplera = " << triples[0][j] << endl;

  //         tmp_gfp_b.input(isb, false);
  //         triples_[1][j]= tmp_gfp_b;
  //         // cout << "P0 triples_b_ = " << tmp_gfp_b << endl;
  //         // cout << "P1 triplerb = " << triples[1][j] << endl;

  //         tmp_gfp_c.input(isc, false);
  //         triples_[2][j]= tmp_gfp_c;
  //         // cout << "P0 triples_c_ = " << tmp_gfp_c << endl;
  //         // cout << "P1 triplerc = " << triples[2][j] << endl;

  //         a_b= (triples[0][j] + triples_[0][j]) * (triples[1][j] + triples_[1][j]);
  //         // cout << "a_b = " << a_b;
  //         c= (triples[2][j] + triples_[2][j]);
  //         // cout << "  c = " << c << endl;
  //         if (a_b != c)
  //           {
  //             cout << "triples error " << endl;
  //             exit(0);
  //           }
  //       }
  //   }
}

void MASCOTTriples::Combine()
{
  /*(1) sampling r, \hat{r}*/
  int currIndex;

  /*(2) combine from multiply to a,c,\hat{a},\hat{c}*/
  // 0-a, 1-b, 2-c, 3-\hat{a}, 4-\hat{c}
  gfp tmp_gfp;
  for (unsigned int j= 0; j < nbTriples; j++)
    {
      combined_triples[0][j].assign_zero(); // a
      combined_triples[2][j].assign_zero(); // c
      combined_triples[3][j].assign_zero(); // a^
      combined_triples[4][j].assign_zero(); // c^

      for (unsigned int jj= 0; jj < tau; jj++)
        {
          currIndex= j * tau + jj;
          tmp_gfp.randomize(G2); // 这个地方怎么保证两方一致的？
          combined_triples[0][j]+= triples[0][currIndex] * tmp_gfp;
          combined_triples[2][j]+= triples[2][currIndex] * tmp_gfp;

          tmp_gfp.randomize(G2);
          combined_triples[3][j]+= triples[0][currIndex] * tmp_gfp;
          combined_triples[4][j]+= triples[2][currIndex] * tmp_gfp;
        }
    }

  // b is already of the right form
  combined_triples[1]= triples[1];
}

void MASCOTTriples::Authenticate()
{

  unsigned int nbPlayers= P.nplayers();
  unsigned int whoami= P.whoami();
  /*Create authenticated shares from the output of Combine*/
  /*Init MAC key and COPE OT*/
  CryptoPP::RandomPool RNG;
  vector<gfp> x_for_COPE(5 * nbTriples + 1);
  gfp x0;
  x0.randomize(G);
  x_for_COPE[0]= x0;
  for (unsigned int i= 0; i < 5 * nbTriples; i++)
    {
      x_for_COPE[i + 1]= combined_triples[(i) % 5][floor((i) / 5)]; // combined_triples里面放的是份额， 拿份额去COPE
    }

  // 诚实方P0 自己生成完整 MAC，并分享出去
  vector<gfp> mac_shares_sender(5 * nbTriples + 1);                 // mac shares 在发送方的份额
  vector<vector<gfp>> q(nbPlayers, vector<gfp>(5 * nbTriples + 1)); // 诚实方会收到，恶意方没用 (a(2) x Delta)(1)
  vector<vector<gfp>> t(nbPlayers, vector<gfp>(5 * nbTriples + 1)); // 恶意方会收到，诚实方没用 (a(2) x Delta)(2)
  if (whoami == 0)                                                  // 只有诚实方P0 拥有完整 MAC key
    {
      if (Delta == 0)
        {
          Delta.randomize(G);
        }
      mac_shares_sender[0]= x_for_COPE[0] * Delta;  // 好像诚实方不需要 x0，但是为了匹配后续接口
      ostringstream os;
      for (unsigned int i= 0; i < 5 * nbTriples; i++)
        {
          gfp complet_mac_shares;
          complet_mac_shares= x_for_COPE[i + 1] * Delta;
          mac_shares_sender[i + 1].randomize(G); // 
          gfp mac_share_receiver= complet_mac_shares - mac_shares_sender[i + 1];

          mac_share_receiver.output(os, false);
          // mac_combined_triples[0][i]+= ;
        }
      P.send_to_player(1, os.str(), 1);
    }
  else
    {
      Delta.assign(0);
    }
  // vector<gfp> mac_shares_sender(5 * nbTriples + 1);
  ostringstream os;
  if (!isInit)
    {
      // for (unsigned int i= 0; i < nbPlayers; i++)
      //   {
      //     if (i < whoami)
      //       {
      //         COPE_R[i].init(P, i, bit_size_field, Delta, RNG, 1);      // COPE_R[0] 1接  0发
      //       }
      //     else if (i > whoami)
      //       {
      //         COPE_S[i].init(P, i, bit_size_field, RNG, 1);          // COPE_S[1]   1发  0接
      //       }
      //   }

      // for (unsigned int i= 0; i < nbPlayers; i++)
      //   {
      //     if (i > whoami)
      //       {
      //         COPE_R[i].init(P, i, bit_size_field, Delta, RNG, 1); // COPE_R[1] 0  1    P0 这边提供Delta
      //       }
      //     else if (i < whoami)
      //       {
      //         COPE_S[i].init(P, i, bit_size_field, RNG, 1);      // COPE_S[0]   1  0
      //       }
      //   }

      // 非对称
      if (whoami == 0)
        {
          COPE_R[1].init(P, 1, bit_size_field, Delta, RNG, 2); // 0 是接收方， 提供Delta
        }
      else
        {
          COPE_S[0].init(P, 0, bit_size_field, RNG, 2); // 1 是发送方，
        }

      cout << "[MASCOTTriples - Authenticate] Finished Base-OTs" << endl;
      isInit= true; // 初始化值初始化一次
    }
  /*Input share to get authenticated share*/
  // 1- Sample x0
  // gfp x0;
  // x0.randomize(G);

  // 2- Create additive sharings - WE ALREADY HAVE A SHARING !  // !!! ???，直接用的已有的份额，不再继续分享了（份额的份额？）

  // 3- Call to COPE.extend
  // vector<vector<gfp>> q(nbPlayers, vector<gfp>(5 * nbTriples + 1));
  // vector<vector<gfp>> t(nbPlayers, vector<gfp>(5 * nbTriples + 1));
  // vector<gfp> x_for_COPE(5 * nbTriples + 1);
  // x_for_COPE[0]= x0;
  // for (unsigned int i= 0; i < 5 * nbTriples; i++)
  //   {
  //     x_for_COPE[i + 1]= combined_triples[(i) % 5][floor((i) / 5)]; // combined_triples里面放的是份额， 拿份额去COPE
  //   }

  // for (unsigned int i= 0; i < nbPlayers; i++)
  //   {
  //     if (i < whoami)
  //       {
  //         COPE_S[i].extend_vec(P, x_for_COPE, t[i]); // 对P1 的x  a(2)进行COPE       COPE_S[0] 1发送  t[0]
  //       }
  //     else if (i > whoami)
  //       {
  //         COPE_R[i].extend_vec(P, q[i]);                                            // COPE_R[1]  0接收  q[1]
  //       }
  //   }

  // for (unsigned int i= 0; i < nbPlayers; i++)
  //   {
  //     if (i > whoami)
  //       {
  //         COPE_S[i].extend_vec(P, x_for_COPE, t[i]); // 这里COPE 不是对 x0...x(5*triple) 都认证吗
  //         // t[i]是COPE 中 P1得到的输出
  //       }
  //     else if (i < whoami)
  //       {
  //         COPE_R[i].extend_vec(P, q[i]); // 从另一个COPE（其他方对自己的输入x进行认证）中得到的q,
  //         // 换个角度，q[i] 是COPE中 p0得到的输出
  //       }
  //   }

  // 非对称
  if (whoami == 0)
    {
      COPE_R[1].extend_vec(P, q[1]); // q[i] 是COPE中 p0得到的输出，向量 所有三元组
    }
  else
    {
      string ss;
      gfp tmp_gfp;
      P.receive_from_player(0, ss, 1);
      istringstream is(ss);
      for (unsigned int i= 0; i < 5 * nbTriples; i++)
        {
          tmp_gfp.input(is, false);
          q[1][i + 1]= tmp_gfp; // 
        }

      COPE_S[0].extend_vec(P, x_for_COPE, t[0]); // 这里COPE 不是对 x0...x(5*triple) 都认证
      // t[i]是COPE 中 P1得到的输出
    }

  // Skip 4, included in COPE
  // 5 - define the MAC shares
  // vector<gfp> mac_shares_sender(5 * nbTriples + 1);
  // for (unsigned int i= 0; i < 5 * nbTriples + 1; i++)
  //   {
  //     if (i == 0)
  //       {
  //         mac_shares_sender[i]= x0 * Delta;    // 为什么 x0 要单独认证?  因为x0不是份额
  //       }
  //     else
  //       {
  //         mac_shares_sender[i]= combined_triples[(i - 1) % 5][floor((i - 1) / 5)] * Delta;   // 本地×
  //       }

  //     for (unsigned int j= 0; j < nbPlayers; j++)
  //       {
  //         if (j != whoami)
  //           {
  //             mac_shares_sender[i]+= t[j][i];  //
  //           }
  //       }
  //   }

  // 非对称 define the MAC shares
  if (whoami == 1)
    {
      for (unsigned int i= 0; i < 5 * nbTriples + 1; i++)
        {
          mac_shares_sender[i]= t[0][i];
        }
    }

  // ******************************************************************************************************************************************
  // 6 - Parties sample random vector for linear combination

  // 7 - Compute y and broadcast it
  // gfp tmp_gfp;
  // vector<gfp> y(nbPlayers); // 线性组合 每个player 一个
  // vector<gfp> m(nbPlayers); // 线性组合对应的MAC   每个player 一个
  // y[whoami].assign(0);

  // for (unsigned int j= 0; j < nbPlayers; j++)
  //   {
  //     m[j].assign(0);
  //     for (unsigned int i= 0; i < 5 * nbTriples + 1; i++)
  //       {
  //         tmp_gfp.randomize(G2);
  //         if (j == whoami) // 我自己这一方才计算y，广播给其他方    例如：whoami == 1
  //           {
  //             if (i == 0)
  //               {
  //                 y[j]+= tmp_gfp * x0;
  //               }
  //             else
  //               {
  //                 y[j]+= tmp_gfp * combined_triples[(i - 1) % 5][floor((i - 1) / 5)];
  //               }

  //             m[j]+= tmp_gfp * mac_shares_sender[i];       //whoami=1,计算出 y[1] m[1],
  //           }
  //         else // 其他各方只计算，我方y对应的mac份额就行
  //           {
  //             m[j]+= tmp_gfp * q[j][i];                //whoami=1 m[0] += tmp_gfp * q[0][i]
  //           }
  //       }
  //   }

  // 非对称
  gfp tmp_gfp;
  vector<gfp> y(nbPlayers); // 线性组合 每个player 一个
  vector<gfp> m(nbPlayers); // 线性组合对应的MAC   每个player 一个
  y[whoami].assign(0);
  m[whoami].assign(0);
  for (unsigned int i= 0; i < 5 * nbTriples + 1; i++)
    {
      tmp_gfp.randomize(G2);
      if (whoami == 1) // 我自己这一方才计算y，广播给其他方
        {
          y[1]+= tmp_gfp * x_for_COPE[i];
          m[1]+= tmp_gfp * mac_shares_sender[i];
        }
      else // 其他各方只计算，y对应的mac份额就行
        {
          m[1]+= tmp_gfp * q[1][i]; // 这里的 m[0]
        }
    }

  // ostringstream os3;
  // y[whoami].output(os3, false);
  // P.send_all(os3.str(), 1); // 广播

  // for (unsigned int i= 0; i < nbPlayers; i++)
  //   {
  //     if (i != whoami)
  //       {
  //         string ss;
  //         P.receive_from_player(i, ss, 1);
  //         istringstream is(ss);
  //         y[i].input(is, false); // 其他方收到广播后记录下来 那一方的y值
  //       }
  //   }

  // 非对称广播
  if (whoami == 1)
    {
      ostringstream os3;
      y[whoami].output(os3, false);
      P.send_all(os3.str(), 1); // 这里不需要广播（先这么跑通再更改）
    }

  if (whoami == 0)
    {
      string ss;
      P.receive_from_player(1, ss, 1);
      istringstream is(ss);
      y[1].input(is, false); // 其他方收到广播后记录下来 那一方的y值
    }

  // cout << "m[0] = " << m[0] << endl;
  // cout << "m[1] = " << m[1] << endl;
  // cout << "y[0] = " << y[0] << endl;
  // cout << "y[1] = " << y[1] << endl;

  // cout << "m[0] - y[0]xDelta = " << m[0] - y[0] * Delta << endl;
  // cout << "m[1] - y[1]xDelta = " << m[1] - y[1] * Delta << endl;

  // Mac check
  gfp check;
  vector<vector<gfp>> sigma(nbPlayers, vector<gfp>(nbPlayers)); // m_i - y * delta_i
                                                                // for (unsigned int i= 0; i < nbPlayers; i++)
                                                                //   {
                                                                //     sigma[i][whoami]= m[i] - y[i] * Delta;
                                                                //     Commit_And_Open(sigma[i], P, true, 1);
                                                                //     check.assign(0);
                                                                //     for (unsigned int j= 0; j < nbPlayers; j++)
                                                                //       {
                                                                //         check+= sigma[i][j];
                                                                //       }
                                                                //     if (check != 0)
                                                                //       {
                                                                //         cout << "ABORT CHECK" << endl;
                                                                //         exit(0);
                                                                //       }
                                                                //   }

  // 非对称 //双方都要check
  for (unsigned int i= 0; i < nbPlayers; i++)
    {
      sigma[i][whoami]= m[i] - y[i] * Delta;
      Commit_And_Open(sigma[i], P, true, 1); // h还要具体看看怎么执行的？
      // cout << "sigma =  " << sigma << endl;

      check.assign(0);
      for (unsigned int j= 0; j < nbPlayers; j++)
        {
          check+= sigma[i][j];
        }
      if (check != 0)
        {
          cout << "ABORT CHECK" << endl;
          exit(0);
        }
    }
  // put shares of mac in mac_combined_triples nicely

  // for (unsigned int i= 0; i < nbTriples; i++)
  //   {
  //     for (unsigned int j= 0; j < 5; j++)
  //       {
  //         mac_combined_triples[j][i]= mac_shares_sender[i * 5 + j + 1]; // mac_shares_sender 中放的本方持有的 （对 本方持有的x份额的）mac份额
  //         // 举例：a(2), mac_shares_sender = a(2) x Delta(2) + (a(2) x Delta(1))(2)
  //         for (unsigned int jj= 0; jj < nbPlayers; jj++)
  //           {
  //             if (jj != whoami)
  //               {
  //                 mac_combined_triples[j][i]+= q[jj][i * 5 + j + 1]; // q[jj] 存放的是(a(1) x Delta(2))(2)
  //                 // a(2) x Delta(2) + (a(2)xDelta(1)(2))(2) + a(1)xDelta(2)(2)
  //                 // + a(1)x Delta(1) + (a(1)xDelta(2))(1) + (a(2)xDelta(1))(1) =

  //                 // 非对称 a(1)  a(1)xDelta -> (a(1)xDelta)(1) + (a(1)xDelta)(2)
  //                 //       a(2)        ->      (a(2)xDelta)(1) + (a(2)xDelta)(2)
  //               }
  //           }
  //       }
  //   }

  for (unsigned int i= 0; i < nbTriples; i++)
    {
      for (unsigned int j= 0; j < 5; j++)
        {
          mac_combined_triples[j][i]= mac_shares_sender[i * 5 + j + 1]; // mac_shares_sender 中放的本方持有的 （对 本方持有的x份额的）mac份额
                                                                        // 举例：a(2), mac_shares_sender = a(2) x Delta(2) + (a(2) x Delta(1))(2)
                                                                        // for (unsigned int jj= 0; jj < nbPlayers; jj++)
                                                                        //   {
                                                                        //     if (jj != whoami)
                                                                        //       {
          mac_combined_triples[j][i]+= q[1][i * 5 + j + 1];             // q[jj] 存放的是(a(1) x Delta(2))(2)
                                                                        // a(2) x Delta(2) + (a(2)xDelta(1)(2))(2) + a(1)xDelta(2)(2)
                                                                        // + a(1)x Delta(1) + (a(1)xDelta(2))(1) + (a(2)xDelta(1))(1) =

          // 非对称 a(1)  a(1)xDelta -> (a(1)xDelta)(1) + (a(1)xDelta)(2)
          //       a(2)        ->      (a(2)xDelta)(1) + (a(2)xDelta)(2)
          //     }
          // }
        }
    }
}

void MASCOTTriples::Sacrifice()
{

  unsigned int whoami= P.whoami();
  unsigned int nbPlayers= P.nplayers();
  // 2 - Linear comb to get rho
  gfp tmp_gfp;
  vector<vector<gfp>> rho(nbTriples, vector<gfp>(nbPlayers));
  vector<gfp> mac_rho(nbTriples);
  vector<vector<gfp>> theta(nbTriples, vector<gfp>(nbPlayers));
  vector<gfp> mac_theta(nbTriples);
  for (unsigned int i= 0; i < nbTriples; i++)
    {

      rho[i][whoami]= (tmp_gfp * combined_triples[0][i]) - combined_triples[3][i];
      mac_rho[i]= (tmp_gfp * mac_combined_triples[0][i]) - mac_combined_triples[3][i];

      theta[i][whoami]= (tmp_gfp * combined_triples[2][i]);
      mac_theta[i]= (tmp_gfp * mac_combined_triples[2][i]);
    }

  // 3 - Open rho
  vector<gfp> rho_opened(nbTriples);
  for (unsigned int i= 0; i < nbTriples; i++)
    {
      ostringstream os;
      vector<string> str_rho_bcast;
      str_rho_bcast.resize(nbPlayers);

      rho[i][whoami].output(os, false);
      str_rho_bcast[whoami]= os.str();

      P.Broadcast_Receive(str_rho_bcast, false, 1);
      rho_opened[i].assign(0);
      for (unsigned int j= 0; j < nbPlayers; j++)
        {
          istringstream is(str_rho_bcast[j]);
          rho[i][j].input(is, false);
          rho_opened[i]+= rho[i][j];
        }
    }

  // 4 - Linear comb to get theta
  for (unsigned int i= 0; i < nbTriples; i++)
    {
      theta[i][whoami]= (theta[i][whoami]) - combined_triples[4][i] - (combined_triples[1][i] * rho_opened[i]);
      mac_theta[i]= (mac_theta[i]) - mac_combined_triples[4][i] - (mac_combined_triples[1][i] * rho_opened[i]);
    }

  // 5 - Check(rho,theta,rho_opened,0)

  for (unsigned int i= 0; i < nbTriples; i++)
    {
      tmp_gfp.randomize(G2);
      theta[i][whoami]= rho_opened[i] * tmp_gfp * Delta;
      theta[i][whoami]= mac_rho[i] * tmp_gfp - theta[i][whoami];
      tmp_gfp.randomize(G2);
      theta[i][whoami]= theta[i][whoami] + mac_theta[i] * tmp_gfp;
    }

  // Mac check y, mac_y
  gfp check;
  for (unsigned int i= 0; i < nbTriples; i++)
    {
      Commit_And_Open(theta[i], P, true, 1);
      // cout << "theta = " << theta << endl;
      check.assign(0);
      for (unsigned int j= 0; j < nbPlayers; j++)
        {
          check+= theta[i][j];
          // cout << " check =   " << check << endl;
        }
      if (check != 0)
        {
          cout << "ABORT CHECK at iteration " << i << endl;
          exit(0);
        }
    }

  // Push this batch of triples into usable ones
  for (unsigned int i= 0; i < nbTriples; i++)
    {
      for (unsigned int j= 0; j < 3; j++)
        {
          usable_triples[j].push_back(combined_triples[j][i]);
          mac_usable_triples[j].push_back(mac_combined_triples[j][i]);
        }
    }

  totalTriplesGenerated+= nbTriples;
}

void MASCOTTriples::Multiplication(vector<gfp> &a, vector<gfp> &mac_a, vector<gfp> &b, vector<gfp> &mac_b, vector<gfp> &out, vector<gfp> &mac_out, vector<vector<list<gfp>>> &u_triples, vector<vector<list<gfp>>> &mac_u_triples, deque<mutex> &mtx)
{

  unsigned int nbMults= a.size();
  vector<vector<gfp>> curr_triples;
  vector<vector<gfp>> mac_curr_triples;
  curr_triples.resize(3);
  mac_curr_triples.resize(3);

  // Get enough triples for all multiplications
  unsigned int nbThreads= u_triples.size();
  unsigned int currsize= 0;
  unsigned int currThread;
  while (nbMults > currsize)
    {
      for (unsigned int i= 0; i < nbThreads; i++)
        {
          currThread= (i + last_thread) % nbThreads;

          mtx[currThread].lock();
          for (unsigned int j= 0; j < 3; j++)
            {

              for (unsigned int k= 0; k < nbTriples; k++)
                {
                  while (u_triples[currThread][j].size() < nbTriples)
                    {
                      mtx[currThread].unlock();
                      usleep(200000);
                      mtx[currThread].lock();
                    }
                  curr_triples[j].push_back((u_triples[currThread][j]).front());
                  mac_curr_triples[j].push_back((mac_u_triples[currThread][j]).front());
                  (u_triples[currThread][j]).pop_front();
                  (mac_u_triples[currThread][j]).pop_front();
                }
            }
          mtx[currThread].unlock();
          currsize= curr_triples[0].size();
          if (currsize >= nbMults)
            {
              last_thread= currThread + 1;
              break;
            }
        }
    }

  // Compute epsilon and rho
  vector<gfp> epsilon;
  vector<gfp> mac_epsilon;
  vector<gfp> rho;
  vector<gfp> mac_rho;

  epsilon.resize(nbMults);
  mac_epsilon.resize(nbMults);
  rho.resize(nbMults);
  mac_rho.resize(nbMults);

  for (unsigned int i= 0; i < nbMults; i++)
    {
      epsilon[i]= a[i] - curr_triples[0][i];
      mac_epsilon[i]= mac_a[i] - mac_curr_triples[0][i];

      rho[i]= b[i] - curr_triples[1][i];
      mac_rho[i]= mac_b[i] - mac_curr_triples[1][i];
    }

  // Open rho and epsilon
  vector<gfp> opened_epsilon;
  vector<gfp> opened_rho;
  gfp tmp;

  opened_epsilon.resize(nbMults);
  opened_rho.resize(nbMults);

  for (unsigned int i= 0; i < nbMults; i++)
    {
      ostringstream os;
      vector<string> str_rho_bcast;

      str_rho_bcast.resize(P.nplayers());
      rho[i].output(os, false);
      str_rho_bcast[P.whoami()]= os.str();
      P.Broadcast_Receive(str_rho_bcast, false, 1);

      opened_rho[i].assign(0);

      for (unsigned int j= 0; j < P.nplayers(); j++)
        {
          istringstream is(str_rho_bcast[j]);
          tmp.input(is, false);
          opened_rho[i]+= tmp;
        }
    }

  for (unsigned int i= 0; i < nbMults; i++)
    {
      ostringstream os;
      vector<string> str_epsilon_bcast;

      str_epsilon_bcast.resize(P.nplayers());
      epsilon[i].output(os, false);
      str_epsilon_bcast[P.whoami()]= os.str();
      P.Broadcast_Receive(str_epsilon_bcast, false, 1);

      opened_epsilon[i].assign(0);

      for (unsigned int j= 0; j < P.nplayers(); j++)
        {
          istringstream is(str_epsilon_bcast[j]);
          tmp.input(is, false);
          opened_epsilon[i]+= tmp;
        }
    }

  // Put epsilon, rho and their macs to list for check later
  for (unsigned int i= 0; i < nbMults; i++)
    {
      opened_values.push_back(opened_rho[i]);
      opened_values.push_back(opened_epsilon[i]);

      mac_opened_values.push_back(mac_rho[i]);
      mac_opened_values.push_back(mac_epsilon[i]);
    }

  // Compute the output
  for (unsigned int i= 0; i < nbMults; i++)
    {
      out[i]= curr_triples[2][i] + opened_epsilon[i] * curr_triples[1][i] + opened_rho[i] * curr_triples[0][i];
      mac_out[i]= mac_curr_triples[2][i] + opened_epsilon[i] * mac_curr_triples[1][i] + opened_rho[i] * mac_curr_triples[0][i] + opened_epsilon[i] * opened_rho[i] * Delta;
      if (P.whoami() == 0)
        {
          out[i]+= opened_rho[i] * opened_epsilon[i];
        }
    }
}

void MASCOTTriples::Check()
{
  unsigned long int nbOpened= opened_values.size();
  cout << "[MASCOTTriples - Check ] Checking " << nbOpened << " opened values " << endl;
  vector<gfp> r;

  // Agree on a random r
  r.resize(nbOpened);
  for (unsigned long int i= 0; i < nbOpened; i++)
    {
      r[i].randomize(G2);
    }

  // compute y and mac_y
  gfp y;
  gfp mac_y;

  y.assign(0);
  mac_y.assign(0);
  for (unsigned long int i= 0; i < nbOpened; i++)
    {
      y+= r[i] * opened_values[i];
      mac_y+= r[i] * mac_opened_values[i];
    }

  // Mac check on y and mac_y
  vector<gfp> sigma;
  gfp check;
  sigma.resize(P.nplayers());

  sigma[P.whoami()]= mac_y - y * Delta;
  Commit_And_Open(sigma, P, true, 1);

  check.assign(0);

  for (unsigned int i= 0; i < P.nplayers(); i++)
    {
      check+= sigma[i];
    }

  if (check != 0)
    {
      cout << "ABORT CHECK OPENED VALUES" << endl;
      throw bad_value();
    }
}

void MASCOTTriples::Open(gfp &a, gfp &mac_a, gfp &out)
{
  ostringstream os;
  vector<string> str_a_bcast;
  gfp tmp;

  str_a_bcast.resize(P.nplayers());
  a.output(os, false);
  str_a_bcast[P.whoami()]= os.str();
  P.Broadcast_Receive(str_a_bcast, false, 1);

  out.assign(0);
  for (unsigned int i= 0; i < P.nplayers(); i++)
    {
      istringstream is(str_a_bcast[i]);
      tmp.input(is, false);
      out+= tmp;
    }

  opened_values.push_back(out);
  mac_opened_values.push_back(mac_a);
}
