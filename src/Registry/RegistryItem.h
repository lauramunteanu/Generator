//____________________________________________________________________________
/*!

\class    genie::RegistryItem

\brief    A templated concrete implementation of the RegistryItemI interface.
          Provides an arbitrary basic type (bool, int, double, string) value
          for RegistryI-type containers.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 06, 2004

*/
//____________________________________________________________________________

#ifndef _REGISTRY_ITEM_H_
#define _REGISTRY_ITEM_H_

#include "RegistryItemI.h"

namespace genie {

template<class T> class RegistryItem;
template<class T>
       ostream & operator << (ostream & stream, const RegistryItem<T> & rec);

template<class T> class RegistryItem : public RegistryItemI {

public:

  RegistryItem() { };
  RegistryItem(T item, bool locked = false);
  RegistryItem(const RegistryItem * ri);
  ~RegistryItem() { }

  const type_info & TypeInfo (void) const { return typeid(fItem); }
  const T &         Data     (void) const { return fItem;         }
  void              Lock     (void)       { fIsLocked = true;     }
  void              UnLock   (void)       { fIsLocked = false;    }
  bool              IsLocked (void) const { return fIsLocked;     }

  void Print(ostream& stream) const;

  friend ostream & operator <<
                        <T>(ostream & stream, const RegistryItem<T> & rec);

private:

  T    fItem;    
  bool fIsLocked;    
};

}      // genie namespace

#endif // _REGISTRY_ITEM_H_
