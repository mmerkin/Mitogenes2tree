Submit-block ::= {
  contact {
    contact {
      name name {
        last "file",
        first "Example",
        middle "",
        initials "",
        suffix "",
        title ""
      },
      affil std {
        affil "Example University",
        div "Department of examples",
        city "Example city",
        country "United Kingdom",
        street "01 Example street",
        email "Example.file@testing.com",
        postal-code "TEST"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "file",
            first "Example",
            middle "",
            initials "",
            suffix "",
            title ""
          }
        }
      },
      affil std {
        affil "Example University",
        div "Department of examples",
        city "Example city",
        country "United Kingdom",
        street "01 Example street",
        postal-code "TEST"
      }
    }
  },
  subtype new
}
Seqdesc ::= pub {
  pub {
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "ALT EMAIL:Example.file@testing.com"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "Submission Title:None"
    }
  }
}
