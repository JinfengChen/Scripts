/////////////////////////////////////////////////////////////////////////
//Copyright (C) 2003 Dekang Lin, lindek@cs.ualberta.ca
//
//Permission to use, copy, modify, and distribute this software for any
//purpose is hereby granted without fee, provided that the above
//copyright notice appear in all copies and that both that copyright
//notice and this permission notice appear in supporting documentation.
//No representations about the suitability of this software for any
//purpose is made. It is provided "as is" without express or implied
//warranty.
//
/////////////////////////////////////////////////////////////////////////

#include "str2idmap.h"
#include "tables.h"

class HmmNode;
class Hmm;

/** A transition between two Hmm nodes. */
struct Transition {
  unsigned long _obs;
  HmmNode* _from;
  HmmNode* _to;

  Transition(HmmNode* from, HmmNode* to, unsigned long obs);
};

/** A node in an Hmm object. */
class HmmNode {
  int _time; // The time slot for this node.
  Hmm* _hmm; // The hmm that this node belongs to
  unsigned long _state; // the state 
  vector<Transition*> _ins; // incoming transitions
  vector<Transition*> _outs;// out going transitions
  double _logAlpha; // alpha_t(s) = P(e_1:t, x_t=s);
  double _logBeta;  // beta_t(s)  = P(e_t+1:T |x_t=s);
  Transition* _psi; // the last transition of the most probable path
                    // that reaches this node
public:
  int time() { return _time; }
  unsigned long state() const { return _state;}
  void logAlpha(double logAlpha) { _logAlpha = logAlpha;}
  double logAlpha() const { return _logAlpha;}
  void logBeta(double logBeta) { _logBeta = logBeta;}
  double logBeta() const { return _logBeta;}
  Transition* psi() { return _psi;}
  void psi(Transition* psi) { _psi = psi;}
  vector<Transition*>& ins() { return _ins;}
  vector<Transition*>& outs() { return _outs;}

  void print();

  HmmNode(int time, unsigned long state, Hmm* hmm) {
    _time = time;
    _state = state;
    _logAlpha = _logBeta = 0;
    _psi = 0;
    _hmm = hmm;
  }
  ~HmmNode();
};

/** The possible states at a particular time slot. */
class TimeSlot : public vector<HmmNode*> {
public:
  ~TimeSlot();
};

/** Pseudo Counts */
class PseudoCounts {
  OneDTable _stateCount;
  TwoDTable _transCount;
  TwoDTable _emitCount;
public:
  OneDTable& stateCount() { return _stateCount;}
  TwoDTable& transCount() { return _transCount;}
  TwoDTable& emitCount()  { return _emitCount;}
  void print(Str2IdMap& str2id);
};

/** An Hmm object implements the Hidden Markov Model. */
class Hmm {
  unsigned long _initState; // the initial state
  TwoDTable _transition;    // transition probabilities
  TwoDTable _emission;      // emission probabilities
  Str2IdMap _str2id;        // mapping between strings and integers
  vector<TimeSlot*> _timeSlots; // the time slots

  double _minLogProb;       // log probabilities lower than this are set to 0

  double getTransProb(Transition* trans);
  double getEmitProb(Transition* trans);

  void forward();  // compute the forward probabilities P(e_1:t, X_t=s)
  void backward(); // compute the backward probabilities P(e_t+1:T | X_t=s)

  /** Re-compute the transition and emission probabilities according
      to the pseudo counts. */
  void updateProbs(PseudoCounts& counts);

  /** Accumulate pseudo counts using the BaumWelch algorithm.  The
      return value is the probability of the observations according to
      the current model.  */
  double getPseudoCounts(PseudoCounts& counts);
public:
  /** Add an observation into the Hmm after the current last time
      slot. The states that have non-zero probability of generating
      the observation will be created. The transition between the new
      states and the states in the previous time slots will also be
      created.*/
  void addObservation(unsigned long obs);

  /** Same as void addObservation(unsigned long obs) above, except
      with a different form of parameter. */
  void addObservation(string obs) { addObservation(_str2id.getId(obs)); }


  /** Read the transition and emission probability tables from the
      files NAME.trans and NAME.emit, where NAME is the value of the
      variable name.*/
  void loadProbs(string name);

  /** Save the transition and emission probability tables into the
      files NAME.trans and NAME.emit, where NAME is the value of the
      variable name. If name is "", both tables are printed on the
      standard output. */
  void saveProbs(string name="");

  /** Read the training data from the input stream. Each line in the
      input stream is an observation sequence. */
  void readSeqs(
    istream& istrm, vector<vector<unsigned long>*>& sequences
    );


  /** Find the state sequence (a path) that has the maximum
      probability given the sequence of observations:
         max_{x_1:T} P(x_1:T | e_1:T);
      The return value is the logarithm of the joint probability 
      of the state sequence and the observation sequence:
        log P(x_1:T, e_1:T)
  */
  double viterbi(vector<Transition*>& path);

  /** return the logarithm of the observation sequence: log P(e_1:T) */
  double obsProb();

  /** Train the model with the given observation sequences using the
      Baum-Welch algorithm. */
  void baumWelch(vector<vector<unsigned long>*>& sequences, int maxIterations);


  /** Conversion between the integer id and string form of states and
      observations. */
  string getStr(unsigned long id) { return _str2id.getStr(id);}
  unsigned long getId(string str) { return _str2id.getId(str);}

  /** Clear all time slots to get ready to deal with another
      sequence. */
  void reset();

  /** Print the states at all time slots and the alpha/beta values at
      these states. */
  void print();

  /** Generate seqs observation sequences according to the model. */
  void genSeqs(ostream& ostrm, int seqs);

  /** Generate an observation sequence with up to maxlen elements
      according to the model. */
  void genSeq(vector<unsigned long>& seq);

  Hmm();
  ~Hmm();
};
